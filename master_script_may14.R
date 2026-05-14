# ══════════════════════════════════════════════════════════════════════════════
#  LIGHT NICHES
#  Predicts underwater PAR at HOBO logger locations
# ══════════════════════════════════════════════════════════════════════════════
#DEM .tif  ──┐
#HOBO .shp ──┼──► prepare_dem / load_hobo_locations / build hourly time series
#env_long  ──┘
#│
#▼
#get_site_centroid()    (lat/lon centre of DEM)
#│
#┌─────────┴────────────────-─┐
#│    for each hour           │
#│                            │
#▼                            ▼
#get_sun()              calc_water_attenuation()
#(altitude + azimuth)       exp(−Kd × depth) → 0–1
#│                            │
#▼                            │
#make_shadow()                │
#ray_shade() → 0–1 matrix     │
#│                            │
#▼                            │
#extract_hobo_shadow()        │
#(sample matrix at points)    │
#│                            │
#└──────────┬──────────────-──┘
#▼
#pred_par
#surface_par × shadow × water_atten
#│
#▼
#Validation
#pred_par ↔ light_mol



# ── LIBRARIES ─────────────────────────────────────────────────────────────────

library(terra)       # raster/vector GIS package
library(raster)      # needed by rayshader
library(rayshader)   # ray_shade()casts shadow maps from a DEM matrix
library(rgl)         # used by rayshader
library(suncalc)     # used for get sunlight position, uses long/lat/time
library(sf)          # read shapefiles, reproject etc.
library(dplyr)       # data wrangling
library(tidyr)       # data reshaping
library(lubridate)   # date-time helper
library(readr)       # tidyverse uses this



# ── WORKING DIRECTORY ─────────────────────────────────────────────────────────
#set working directory, need to change to github repository 

setwd("~/Desktop/MEM MSc/1thesis/LightNiches-main/output")


# ── LOAD INPUT DATA ───────────────────────────────────────────────────────────

#env_long datasheet, and converts timestamps into a date-time object in the correct timezone for lizard island
#this sheet has kd back-calculated from the sensors, the code doesnt use it, I need to remove this
env_long <- read.csv("env_long_depth_kd_longlat.csv") %>%
  mutate(datetime_pos = as.POSIXct(datetime_pos, tz = "Australia/Brisbane"))


# geom_new.csv provides correct lon/lat (in local transmercator projection) and
# new_unit (= hoboID) for every logger at every site/year.
# this replaces shapefile names, which didnt have names that aligned and were never correctly joined to env_long.
geom <- read.csv("geom_new.csv") %>%
  distinct(rec, year, new_unit, lon, lat) %>%     # one row per logger per site/year
  rename(site = rec, hoboID = new_unit) %>% #renames columns to match with env_long
  mutate(year = as.integer(year))


# ── FILE-PATH ROOTS ───────────────────────────────────────────────────────────
#create filepath for DEMS 
DEM_DIR <- "/Users/samueleszenyi/Desktop/MEM MSc/1thesis/dems"

# ── Kd — diffuse attenuation coefficient (m⁻¹)  ← change this number ────────
#stores the universal KD number, I need to adjust this a bit better, 
#look at jervlov water quality values etc. Right now (May 14) I am using a mean value from a site

KD_CONSTANT <- 0.31


# ══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# reusable code chunk recipes that I made to be called later 
# ══════════════════════════════════════════════════════════════════════════════

# ── trim_na_border ────────────────────────────────────────────────────────────
# Remove rows/columns whose NA fraction exceeds threshold.
trim_na_border <- function(mat, threshold = 0.5) { #takes a grid matrix and removes rows or columns where more than half the values are missing 
  #I can change the threshold later 
  keep_rows <- apply(mat, 1, function(r) mean(is.na(r)) <= threshold)
  mat <- mat[keep_rows, , drop = FALSE]
  keep_cols <- apply(mat, 2, function(c) mean(is.na(c)) <= threshold)
  mat <- mat[, keep_cols, drop = FALSE]
  mat
}

# ── prepare_dem ───────────────────────────────────────────────────────────────
# Load DEM, downsample, convert to a matrix for rayshading.
prepare_dem <- function(site, year, fact = 3) { #the downsample is by a factor of 3, and makes raytraycing faster without losing a lot of noise
  path <- file.path(DEM_DIR, paste0(site, year, ".tif"))
  message("  Loading DEM: ", path)
  dem        <- rast(path) |> trim()
  dem_small  <- aggregate(dem, fact = fact)
  # zscale = 1: horizontal and vertical units are both metres, so there is no crazy exaggeerations
  zscale     <- 1
  dem_matrix <- raster_to_matrix(raster(dem_small))
  dem_matrix <- trim_na_border(dem_matrix, threshold = 0.5)
  message("  DEM loaded: ", nrow(dem_matrix), " x ", ncol(dem_matrix),
          " | raster: ", nrow(dem_small), " x ", ncol(dem_small))
  list(matrix = dem_matrix, rast = dem_small, zscale = zscale)
}
#this returns 3 things, a matrix for rayshader, a raster for coordinate extraction, and the zscale

# ── load_hobo_locations ───────────────────────────────────────────────────────
# build a sf point object from geom_new coordinates for a given site/year.
# coords are in the local transverse mercator projection used when the
# geom script vivi made extracted them from the DEMs.

load_hobo_locations <- function(site, year, geom) {
  tm_crs <- paste0("+proj=tmerc +lat_0=-14.6989736557007", #this defines the coordinate system the numbers are in 
                   " +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0",
                   " +datum=WGS84 +units=m +no_defs")
  locs <- geom %>%
    filter(site == !!site, year == !!year) %>%
    filter(!is.na(lon), !is.na(lat)) %>%
    distinct(hoboID, .keep_all = TRUE)

  if (nrow(locs) == 0) {
    stop("No logger locations found in geom for site=", site, " year=", year)
  }

  message("  Loaded ", nrow(locs), " logger locations from geom for ",
          site, year)
  st_as_sf(locs, coords = c("lon", "lat"), crs = tm_crs)
}

# ── get_sun ───────────────────────────────────────────────────────────────────
#this uses a timestamp, lat and long to get the sun altitude (height above horizon in radians)
#and azimuth, which is the compass direction in radians
get_sun <- function(time, lat, lon) {
  sun <- getSunlightPosition(date = time, lat = lat, lon = lon)
  list(altitude = sun$altitude, azimuth = sun$azimuth)
}
#this determines the angle / strength of the shadows

# ── make_shadow ───────────────────────────────────────────────────────────────
#this runs the raytracing in the background, sending the dem_matrix that I made with prepare_dem
#and sun positioning to do raytracing using rayshader, this returns a grid of values 
#0 is fully shaded, and 1 is fully lit, if the sun is below the horizon it skips the measurements to save time
make_shadow <- function(dem_matrix, sun, zscale) {
  if (sun$altitude <= 0.1) {
    message("    (sun below horizon — skipping raytrace)")
    return(NULL)
  }
  shadow <- ray_shade(
    dem_matrix,
    sunaltitude = sun$altitude * 180 / pi,
    sunangle    = sun$azimuth  * 180 / pi,
    zscale      = zscale
  )
  message("    shadow range: ", round(min(shadow, na.rm = TRUE), 3),
          " to ", round(max(shadow, na.rm = TRUE), 3),
          " | mean: ", round(mean(shadow, na.rm = TRUE), 3))
  shadow
}

# ── calc_water_attenuation ────────────────────────────────────────────────────
#Beer-Lambert law: exp(−Kd × depth)
#this applies beer lambert law, returning a fraction between 0-1 representing how much light survives the trip through the water column to a certain point
#a deeper logger would get a smaller fraction
calc_water_attenuation <- function(Kd, water_depth_m) {
  ifelse(
    is.na(Kd) | water_depth_m <= 0,
    NA_real_,
    exp(-Kd * water_depth_m)
  )
}

# ── extract_hobo_dem ──────────────────────────────────────────────────────────
#this looks up hobo logger locations and stores it as a new column called dem elevation, which is not used in the final prediction but good to have 
extract_hobo_dem <- function(dem_rast, hobo_sf) {
  hobo_proj        <- st_transform(hobo_sf, crs(dem_rast))
  vals             <- terra::extract(dem_rast, vect(hobo_proj))
  hobo_sf$dem_elev <- vals[, 2]
  hobo_sf
}

# ── extract_hobo_shadow ───────────────────────────────────────────────────────
#this takes the shadow grid from make_shadow, and puts it in to a spatial raster with the same coordinates as the raster, than samples each loggers exact location
#to get the logger's shadow fraction for that hour
extract_hobo_shadow <- function(shadow_matrix, dem_rast, hobo_sf) {
  shadow_rast      <- rast(shadow_matrix)
  ext(shadow_rast) <- ext(dem_rast)
  crs(shadow_rast) <- crs(dem_rast)
  hobo_proj        <- st_transform(hobo_sf, crs(dem_rast))
  vals             <- terra::extract(shadow_rast, vect(hobo_proj))
  vals[, 2]
}

# ── get_site_centroid ─────────────────────────────────────────────────────────
# this finds the center point of the DEM in regular lat long degrees (WGS 84/EPSG 4326)
#this is used for the solar position in the main loops
get_site_centroid <- function(dem_rast) {
  e <- ext(project(dem_rast, "EPSG:4326"))
  list(lat = (e$ymin + e$ymax) / 2, lon = (e$xmin + e$xmax) / 2)
}


# ══════════════════════════════════════════════════════════════════════════════
#  PER-TIMESTEP PROCESSING
#the core calculations for each timestep (I chose 1 hour)
# ══════════════════════════════════════════════════════════════════════════════
process_timestep <- function(env_block, dem_matrix, dem_rast, dem_zscale,
                             hobo_sf, centroid) {

  # 1. calculate solar position for the chosen hour at the chosen site
  timestamp <- env_block$datetime_pos[1]
  sun <- get_sun(timestamp, centroid$lat, centroid$lon)

  message("  >> Raytracing: ", format(timestamp, "%Y-%m-%d %H:%M"),
          " | sun alt: ", round(sun$altitude * 180 / pi, 1), "\u00b0",
          " | sun az: ",  round(sun$azimuth  * 180 / pi, 1), "\u00b0",
          " | Kd: ",      round(env_block$Kd[1], 4))

  # 2 get the topographic shading from the dem
  shadow <- make_shadow(dem_matrix, sun, zscale = dem_zscale)

  # 3 gets shadow values for each hobo logger, if the shadow map exists, it samples each logger location
  #if the sun is below the horizon (dark hours) i give every logger an NA and skip 
  shadow_vals <- if (!is.null(shadow)) {
    extract_hobo_shadow(shadow, dem_rast, hobo_sf)
  } else {
    rep(NA_real_, nrow(hobo_sf))
  }

  # 4. this makes a table pairing each shadow value with the logger id
  shadow_df <- data.frame(
    hoboID_env      = hobo_sf$hoboID,
    shadow_fraction = shadow_vals
  )

  # 5 this joins the shadow fraction with the loger data (depth, surface par, kd) and then computes the final predicted value as 
  # surface par x shadow fraction x water attenuation 
  result <- shadow_df %>%
    left_join(
      env_block %>%
        mutate(hoboID_env = as.character(hoboID)) %>%
        select(hoboID_env, Kd, water_depth_m, surface_par, site, year),
      by = "hoboID_env"
    ) %>%
    mutate(
      datetime_pos = timestamp,
      sun_altitude = sun$altitude,
      sun_azimuth  = sun$azimuth,
      water_atten  = calc_water_attenuation(Kd, water_depth_m),
      pred_par     = surface_par * shadow_fraction * water_atten
    )

  as_tibble(result)
}


# ══════════════════════════════════════════════════════════════════════════════
#  main code loop that puts it all together
# master function that does following:

# loops over every site x year combination in data
# loads elevation map and logger positions ONCE
#builds a shcedule of hourly daytime timesteps to process
      #and for every hour that qualifies (sun over horizon), raytraces DEM and computes PAR at every logger location

#the outer loop sets up the machinery for every site/year, and the inner loop runs it one hour at a time

#inputs, env_long - full environmental data table 
#geom - logger coorinate lookup table
#dems

#output - one table with one row per logger per timestep, with pred_par
# ══════════════════════════════════════════════════════════════════════════════

run_pipeline <- function(env_long, geom) {
  #results in an empty list that it fills one bit at a time

  results <- list()
  #get the unique site codes present (CB, L1, Gt, etc)
  sites   <- unique(env_long$site)
  #unique skips duplicates so we only do one site once 

  message("=================================================")
  message("Starting pipeline")
  message("Sites to process: ", paste(sites, collapse = ", "))
  message("=================================================")

#outer loop 1 -----------------------------------------------------------------
  #each site ha its own dem file and logger positions, so it processes them one at a time
  for (s in sites) {
    years <- unique(env_long$year[env_long$site == s])
#find which years of data exist for each site
#the env_long$year[..... ] means give me the year column but for only rows where the site column equals s
    
#outer loop 2 -----------------------------------------------------------------
#---------step1 load spatial data 
    #iterate over years within site:
    #load spatial data once per site/year
    #when they are outside hourly loop we only do this once, since it takes a lot of time (i broke my R doing it inside the loop)
    for (y in years) {
      message("\n--- Site: ", s, " | Year: ", y, " ---")
#load dem file, downsample it, convert to matrix
      #returns $matrix, $rast, $zscale
      # ── Load spatial data once per site/year ──────────────────────────────
      dem      <- prepare_dem(s, y) #build spatial object of logger locations with HOBO ID that matches env_long 
      hobo_sf  <- load_hobo_locations(s, y, geom)   # uses geom_new, not shapefilesto sample DEM elevation at logger locations
      hobo_sf  <- extract_hobo_dem(dem$rast, hobo_sf) 
      centroid <- get_site_centroid(dem$rast) #find geographic center of DEM to use for sunlight, the sites are small enough to use once center
      message("  Site centroid: lat ", round(centroid$lat, 4),
              " lon ", round(centroid$lon, 4))
#---------step2 build hourly schedule of timesteps
      #the raw data is 5 minute intervals, but we only want to raytrace 1 hour or 30 minutes to save computing
      #the section here collapses 5 min data to one row per loggere per hour, computes sun angle, and throws away hours where the sun is too low to mean anything
      #%>% chains everything together in order 
      hourly_times <- env_long %>%
        
        #keep only rows belonging to this site year
        filter(site == s, year == y) %>%
        
        #sort chronologically for floor date to work
        arrange(datetime_pos) %>%
        
        #snap every 5 min timestamp to its hour 
        #11:05, 11:10, 11:55 all become 11, which makes rows share the same datetime
        mutate(datetime_pos = floor_date(datetime_pos, "1 hour")) %>%
        
        #group rows by hgour and logger id, so each hour/logger gets its own independent group
        group_by(datetime_pos, hoboID) %>%
        
        #keep only the first row, since we just want one representative row per logger per hour
        slice(1) %>%
        ungroup() %>%
        
        #drop any remaining duplicate rows and keep only columns that process_timestep needs 
        distinct(datetime_pos, hoboID, site, year,
                 surface_par, tide_cm, water_depth_m) %>%
        
        #give every row with our constant KD value, the kd values from env_long are discarded, since using them would be circular
        mutate(Kd = KD_CONSTANT) %>%     
        
        #r applies the next mutate one row at a time, getsunlightposition needs one single timestapm
        rowwise() %>%
        
        #computer solar altitue in degrees, stores it as sun_alt, conveting radians to degrees by altitude*180/pi 
        mutate(
          sun_alt = getSunlightPosition(
            datetime_pos, centroid$lat, centroid$lon
          )$altitude * 180 / pi #return degrees
        ) %>%
        ungroup() %>%
        
        #discard hours where the sun is less than 15 degrees over the horizon, this can also be changed
        #at supr low angles, the shadows hit dem mperfections harder and the raytracer produces unreliable shadows
        #we lose a bit of early morning data but it can always be changed
        filter(sun_alt > 15) %>%
        
        #remove sun_alt column, since it served its purpose
        select(-sun_alt)

      #get list of hours to raytrace, we only need one shadow map per hours
      timestamps <- unique(hourly_times$datetime_pos)
      
      #count how many steps we calculate, which is ussed for the message below
      n_steps    <- length(timestamps)

      message("  Timesteps to process: ", n_steps,
              " (", format(min(timestamps), "%Y-%m-%d %H:%M"),
              " to ", format(max(timestamps), "%Y-%m-%d %H:%M"), ")")

#inner loop -----------------------------------------------------------------
#---------step 3 inner loop, one raytracing per hour
      #for each hour i pull out all logger rows for that hour, pass them to process_timestep, which raytraces once than extracts shaodwvalues, then output results
      
#forseq_len(n_steps) produces 1,2,3,n_steps, which is safer than writing 1:n_steps because if n_steps were 0, 1:0 would loop twice
      for (i in seq_len(n_steps)) {
        
        #pick the i-th timestamp from the hours we decided
        t         <- timestamps[i]
        
        #pull out logger rows for this hour, [] is a row filter, which keeps rows where datetime_pos == t
        #processtimestep needs all loggers at once so it can extract shadow values for each one from the shaodw map that was made earlier
        env_block <- hourly_times[hourly_times$datetime_pos == t, ]

        message("\n[", i, "/", n_steps, "] site: ", s, " | year: ", y,
                " | loggers in block: ", nrow(env_block),
                " | Kd: ", KD_CONSTANT)

        #run the full prediction for the timestep
        #solar position -> raytrace -> logger shadows -> beer lamber law -> predicted par (pred_par)
        res <- process_timestep(
          env_block, #data frame: one row per logger for this hour
          dem$matrix,#dem as a numeric matrix
          dem$rast, #dem as a spatial object/raster for coordinate stuff
          dem$zscale, #z-scale = 1 (m horizontal = m vertical, so there is no crazy exaggerations)
          hobo_sf, #logger locations as a spatial object
          centroid #site lat/lon for suncalc sun angle lookup 
        )

        #add hour's result tibble to results
        #length(results) +1 puts it as the next empty row 
        results[[length(results) + 1]] <- res
      }

      message("\n  Done with site: ", s, " | year: ", y)
    }
  }

  message("\n=================================================")
  message("Pipeline complete :) Total result chunks: ", length(results))
  message("=================================================\n")

  #combine all hour's tibbles into one single df and return it
  #creating the final output of the data pipeline
  bind_rows(results)
}


# ══════════════════════════════════════════════════════════════════════════════
#  RUN a test/chunk of it— edit site/year filter here
# ══════════════════════════════════════════════════════════════════════════════

env_long_test <- env_long %>%
  filter(site == "L1", year == 19)

output <- run_pipeline(env_long_test, geom)


# ══════════════════════════════════════════════════════════════════════════════
#  VALIDATION — join predicted PAR to measured PAR
# ══════════════════════════════════════════════════════════════════════════════

validation <- output %>%
  left_join(
    env_long_test %>%
      mutate(hoboID_env = as.character(hoboID)) %>%
      select(datetime_pos, hoboID_env, light_mol),
    by = c("datetime_pos", "hoboID_env")
  ) %>%
  filter(!is.na(pred_par), !is.na(light_mol))

nrow(validation)

write.csv(validation, "L119_output.csv", row.names = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT — measured vs predicted PAR per logger
# ══════════════════════════════════════════════════════════════════════════════

hobos   <- sort(unique(validation$hoboID_env))
n_hobos <- length(hobos)
n_cols  <- 5
n_rows  <- ceiling(n_hobos / n_cols) + 1   # +1 row for legend

par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 1.5, 0.5))

for (id in hobos) {
  sub  <- validation[validation$hoboID_env == id, ]
  sub  <- sub[order(sub$datetime_pos), ]
  ylim <- range(c(sub$pred_par, sub$light_mol), na.rm = TRUE)

  plot(sub$datetime_pos, sub$light_mol,
       type = "l", col = "blue", lwd = 1.5,
       ylim = ylim, xlab = "", ylab = "",
       main = id, cex.main = 0.8, xaxt = "n")

  axis.POSIXct(1, sub$datetime_pos, format = "%H", cex.axis = 0.6)

  lines(sub$datetime_pos, sub$pred_par,
        col = "red", lwd = 1.5, lty = 2)
}

plot.new()
legend("center",
       legend = c("Measured", "Predicted"),
       col    = c("blue", "red"),
       lty    = c(1, 2), lwd = 2,
       cex    = 0.9, bty = "n")

par(mfrow = c(1, 1))

env_long_test %>%
  mutate(hour = lubridate::floor_date(datetime_pos, "1 hour")) %>%
  group_by(hour) %>%
  summarise(surface_par = first(surface_par)) %>%
  arrange(desc(surface_par)) %>%
  head(10)


# 1. How many rows did the pipeline actually produce?
nrow(output)

# 2. What hoboIDs are in output vs env_long_test?
sort(unique(output$hoboID_env))
sort(unique(as.character(env_long_test$hoboID)))

# 3. Did the join match anything before the filter?
output %>%
  left_join(
    env_long_test %>%
      mutate(hoboID_env = as.character(hoboID)) %>%
      select(datetime_pos, hoboID_env, light_mol),
    by = c("datetime_pos", "hoboID_env")
  ) %>%
  summarise(
    total = n(),
    has_pred = sum(!is.na(pred_par)),
    has_light = sum(!is.na(light_mol)),
    has_both = sum(!is.na(pred_par) & !is.na(light_mol))
  )





