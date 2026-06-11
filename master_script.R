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
#           ▼
#       pred_par
#surface_par × shadow × water_atten
#   │
#   ▼
# Validation
# pred_par ↔ light_mol



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
  mutate(datetime_pos = as.POSIXct(datetime_pos, format = "%m/%d/%y %H:%M", tz = "Australia/Brisbane"))


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

KD_CONSTANT <- 0.15
SHADOW_GAMMA <- 1
Z_SCALE <- 0.006
TIMESTEP_INTERVAL <- "1 hour"
DOWNSAMPLE_FACTOR <- 3
MAXSEARCH_FACT <- 400

# ══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# reusable code chunk recipes that I made to be called later 
# ══════════════════════════════════════════════════════════════════════════════

# ── trim_na_border ────────────────────────────────────────────────────────────
# Remove rows/columns whose NA fraction exceeds threshold.
trim_na_border <- function(mat, threshold = 0.75) { #takes a grid matrix and removes rows or columns where more than half the values are missing 
  #I can change the threshold later 
  keep_rows <- apply(mat, 1, function(r) mean(is.na(r)) <= threshold)
  mat <- mat[keep_rows, , drop = FALSE]
  keep_cols <- apply(mat, 2, function(c) mean(is.na(c)) <= threshold)
  mat <- mat[, keep_cols, drop = FALSE]
  mat
}

# ── prepare_dem ───────────────────────────────────────────────────────────────
# Load DEM, downsample, convert to a matrix for rayshading.
prepare_dem <- function(site, year, fact = DOWNSAMPLE_FACTOR) {
  path <- file.path(DEM_DIR, paste0(site, year, ".tif"))
  message("  Loading DEM: ", path)
  dem        <- rast(path) |> trim()
  message("  Raw DEM CRS: ", crs(dem, describe=TRUE)$name)
  message("  Raw DEM extent: ", paste(as.vector(ext(dem)), collapse = " "))
  dem_small  <- aggregate(dem, fact = fact)
  zscale     <- Z_SCALE
  dem_matrix <- raster_to_matrix(raster(dem_small))
  
  # trim the matrix AND crop the raster to match
  keep_rows <- apply(dem_matrix, 1, function(r) mean(is.na(r)) <= 0.5)
  keep_cols <- apply(dem_matrix, 2, function(c) mean(is.na(c)) <= 0.5)
  dem_matrix <- dem_matrix[keep_rows, keep_cols, drop = FALSE]
  dem_small  <- dem_small[keep_rows, keep_cols, drop = FALSE]
  
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

#slow command
make_shadow <- function(dem_matrix, sun, zscale, gamma = SHADOW_GAMMA) { # gamma raises every shadow fraction (0–1) to a power of (gamma = x), which creates a curve where partly-shaded/shaded areas are pushed darker while fully-lit areas near or stay at 1
  if (sun$altitude <= 0.1) {
    message("    (sun below horizon — skipping raytrace)")
    return(NULL)
  }
  shadow  <- ray_shade(
    dem_matrix,
    sunaltitude = sun$altitude * 180 / pi, #gets sun altitude in degrees from radians
    sunangle    = sun$azimuth  * 180 / pi, #gets sun azimuth in degrees from radians
    zscale      = Z_SCALE, #the zscale remains the same (1:1 horizontal vs vertical)
    maxsearch   = MAXSEARCH_FACT #limits how each ray travels, so smaller features have more importance (but this takes longer)
  )
  ambient <- ambient_shade(dem_matrix, zscale = zscale) #calculates a new shadow layer based on ambient shade, making more enclosed spots darker regardless of sun angle
  shadow  <- (shadow * ambient) ^ gamma #darkens areas by applying a power of gamma to the shadow*ambient layer,this will crush the shadow layer (hopefully) to make darker stronger and lighter lighter 
  message("    shadow range: ", round(min(shadow, na.rm = TRUE), 3),
          " to ", round(max(shadow, na.rm = TRUE), 3),
          " | mean: ", round(mean(shadow, na.rm = TRUE), 3))
  shadow
}

#original command 
make_shadow <- function(dem_matrix, sun, zscale) {
  if (sun$altitude <= 0.1) {
    message("    (sun below horizon — skipping raytrace)")
    return(NULL)
  }
  shadow <- ray_shade(
    dem_matrix,
    sunaltitude = sun$altitude * 180 / pi,
    sunangle    = sun$azimuth  * 180 / pi,
    zscale      = zscale,
    maxsearch = MAXSEARCH_FACT
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
  hobo_wgs84       <- st_transform(hobo_sf, "EPSG:4326")
  vals             <- terra::extract(dem_rast, vect(hobo_wgs84))
  hobo_sf$dem_elev <- vals[, 2]
  hobo_sf
}
# ── extract_hobo_shadow ───────────────────────────────────────────────────────
#this takes the shadow grid from make_shadow, and puts it in to a spatial raster with the same coordinates as the raster, than samples each loggers exact location
#to get the logger's shadow fraction for that hour
extract_hobo_shadow <- function(shadow_matrix, dem_rast, hobo_sf) {
  shadow_rast      <- rast(shadow_matrix)
  ext(shadow_rast) <- ext(dem_rast)
  crs(shadow_rast) <- crs(dem_rast)   # WGS84
  
  # reproject hobo points TO WGS84 to match the shadow raster
  hobo_wgs84 <- st_transform(hobo_sf, "EPSG:4326")
  
  # buffer 10cm = 0.10m, convert to degrees at this latitude
  # at -14.7 degrees: 1 degree ≈ 111,000m, so 0.4m ≈ 0.0000036 degrees
  buffer_deg <- 0.10 / 111000
  hobo_buffered <- st_buffer(hobo_wgs84, dist = buffer_deg)
  
  # extract mean shadow value within each buffer
  vals <- terra::extract(shadow_rast, vect(hobo_buffered), fun = mean, na.rm = TRUE)
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

  message("    shadow per logger: ",
          paste(hobo_sf$hoboID, "=", round(shadow_vals, 3), collapse = ", "))
  message("    dem_rast extent: ", paste(round(as.vector(ext(dem_rast)), 1), collapse = " "))
  hobo_proj <- st_transform(hobo_sf, crs(dem_rast))
  message("    hobo bbox: ", paste(round(st_bbox(hobo_proj), 1), collapse = " "))
  message("    hobo coords raw: ", paste(round(st_coordinates(hobo_sf), 1), collapse = " "))
  message("    dem_rast CRS: ", crs(dem_rast, describe=TRUE)$name)
  message("    hobo_sf CRS: ", st_crs(hobo_sf)$input)
  
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
      pred_par     = surface_par * (shadow_fraction ^ SHADOW_GAMMA) * water_atten
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
        mutate(datetime_pos = floor_date(datetime_pos, TIMESTEP_INTERVAL)) %>%
        
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
#env_long_hour <- env_long %>%
#  filter(site == "CB", year == 18) %>%
#  filter(datetime_pos >= as.POSIXct("2018-12-04 11:00:00", tz = "Australia/Brisbane"),
#         datetime_pos <= as.POSIXct("2018-12-04 13:00:00", tz = "Australia/Brisbane"))


env_long_test <- env_long %>%
  filter(site == "CB", year == 18)

env_long_test <- env_long %>%
  filter(site == "CB", year == 18) %>%
  filter(datetime_pos >= as.POSIXct("2018-12-04 11:00:00", tz = "Australia/Brisbane"),
         datetime_pos <  as.POSIXct("2018-12-04 18:00:00", tz = "Australia/Brisbane"))

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

#write.csv(validation, "L119_output.csv", row.names = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT — measured vs predicted PAR per logger
# ══════════════════════════════════════════════════════════════════════════════
# ── PLOT ──────────────────────────────────────────────────────────────────────
dev.off()
dev.new()

site_name  <- unique(validation$site)
date_range <- format(range(validation$datetime_pos), "%Y-%m-%d")
plot_title <- paste0(
  site_name, " | ", date_range[1], " to ", date_range[2], "\n",
  "Kd=", KD_CONSTANT,
  " | gamma=", SHADOW_GAMMA,
  " | zscale=", Z_SCALE,
  " | interval=", TIMESTEP_INTERVAL
)

ylim_shared <- range(c(validation$pred_par, validation$light_mol), na.rm = TRUE)
hobos       <- sort(unique(validation$hoboID_env))
n_hobos     <- length(hobos)
n_cols      <- 5
n_rows      <- ceiling(n_hobos / n_cols)
legend_cell <- n_hobos + 2

# pad with 0 (= empty cell) instead of NA
plot_cells <- c(2:(n_hobos + 1), rep(0, n_cols * n_rows - n_hobos))

layout_mat <- rbind(
  rep(1, n_cols),
  matrix(plot_cells, nrow = n_rows, ncol = n_cols, byrow = TRUE),
  rep(legend_cell, n_cols)
)

layout(layout_mat, heights = c(0.8, rep(3, n_rows), 0.8))

# 1. title
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, plot_title, cex = 1.2, font = 2, adj = 0.5)

# 2. per-logger panels
for (i in seq_along(hobos)) {
  par(mar = c(2, 2, 1.5, 0.5))
  id  <- hobos[i]
  sub <- validation[validation$hoboID_env == id, ]
  sub <- sub[order(sub$datetime_pos), ]
  
  is_left <- (i - 1) %% n_cols == 0
  
  plot(sub$datetime_pos, sub$light_mol,
       type  = "l", col = "blue", lwd = 1.5,
       ylim  = ylim_shared,
       xlab  = "", ylab = "",
       main  = id, cex.main = 0.85,
       xaxt  = "n",
       yaxt  = if (is_left) "s" else "n")
  
  axis.POSIXct(1, sub$datetime_pos, format = "%H", cex.axis = 0.6)
  lines(sub$datetime_pos, sub$pred_par, col = "red", lwd = 1.5, lty = 2)
}

# 3. legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("Measured", "Predicted"),
       col    = c("blue", "red"),
       lty    = c(1, 2), lwd = 2,
       cex    = 1.0, bty = "n",
       horiz  = TRUE)

layout(1)

#Shadow map Hourly (simple) with HOBOS on it to see if they are in correct places, zscale exaggerated
dem <- prepare_dem("CB", 18)
hobo_sf <- load_hobo_locations("CB", 18, geom)
hobo_wgs84 <- st_transform(hobo_sf, "EPSG:4326")

shadow_test <- ray_shade(dem$matrix,
                         sunaltitude = 72,
                         sunangle    = -62.8,
                         zscale      = 0.005,
                         maxsearch   = 200)

shadow_rast <- rast(shadow_test)
ext(shadow_rast) <- ext(dem$rast)
crs(shadow_rast) <- crs(dem$rast)

#png("shadow_overlay.png", width = 1200, height = 1200)
plot(shadow_rast, main = "Shadow map with logger locations")
plot(vect(hobo_wgs84), add = TRUE, col = "red", cex = 2, pch = 16)
text(st_coordinates(hobo_wgs84)[,1], 
     st_coordinates(hobo_wgs84)[,2], 
     labels = hobo_sf$hoboID, 
     col = "white", cex = 0.7, font = 2,
     pos = 3)  # pos=3 puts label above the point
dev.off()



#plot of DLI with sensors
# calculate DLI from measured data (5 min intervals)
dli_measured <- aggregate(
  light_mol ~ hoboID + as.Date(datetime_pos),
  data = env_long[env_long$site == "CB" & env_long$year == 18, ],
  FUN = function(x) sum(x, na.rm=TRUE) * (5/60) / 1e6
)
colnames(dli_measured) <- c("hoboID", "date", "dli_measured")

# calculate DLI from predicted (hourly)
validation$date <- as.Date(validation$datetime_pos)
dli_predicted <- aggregate(
  pred_par ~ hoboID_env + date,
  data = validation,
  FUN = function(x) sum(x, na.rm=TRUE) / 1e6
)
colnames(dli_predicted) <- c("hoboID", "date", "dli_predicted")

# join
dli <- merge(dli_measured, dli_predicted, by = c("hoboID", "date"))

# plot
plot(dli$dli_measured, dli$dli_predicted,
     pch = 16, col = "steelblue",
     xlab = "measured DLI (mol m-2 d-1)",
     ylab = "predicted dli (mol m-2 d-1)",
     main = "daly light integral, measured vs predicted, CB18")

abline(0, 1, lty = 2, col = "red", lwd = 2)

text(dli$dli_measured, dli$dli_predicted,
     labels = dli$hoboID,
     cex = 0.7, pos = 3)


###### Plot of just points at different times of day 

# filter to specific hours
time_df <- validation %>%
  mutate(hour = as.integer(format(datetime_pos, "%H"))) %>%
  filter(hour %in% c(12, 15, 17)) %>%
  mutate(time_label = case_when(
    hour == 12 ~ "12:00",
    hour == 15 ~ "15:00",
    hour == 17 ~ "17:00"
  ))

# get shared axis limits
lim <- range(c(time_df$light_mol, time_df$pred_par), na.rm = TRUE)

# plot
plot(NULL, xlim = lim, ylim = lim,
     xlab = "Observed PAR (µmol)",
     ylab = "Predicted PAR (µmol)",
     main = "Predicted vs Observed PAR by time of day, CB18, KD=0.2, zscale exaggerated")

abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)

# add points per time with different symbols
points(time_df$light_mol[time_df$hour == 12],
       time_df$pred_par[time_df$hour == 12],
       pch = 16, col = "orange", cex = 1.5)

points(time_df$light_mol[time_df$hour == 15],
       time_df$pred_par[time_df$hour == 15],
       pch = 15, col = "steelblue", cex = 1.5)

points(time_df$light_mol[time_df$hour == 17],
       time_df$pred_par[time_df$hour == 17],
       pch = 8, col = "darkred", cex = 1.5)

legend("topleft",
       legend = c("12:00", "15:00", "17:00"),
       pch = c(16, 15, 8),
       col = c("orange", "steelblue", "darkred"),
       bty = "n")


# load and plot the DEM with a terrain colour palette
dem_view <- rast("/Users/samueleszenyi/Desktop/MEM MSc/1thesis/dems/CB18.tif")

plot(dem_view, 
     col = terrain.colors(100),
     main = "CB18 DEM")



######DEM ZSCALE CALCULATION

# calculate correct zscale from DEM properties
dem_check <- rast("/Users/samueleszenyi/Desktop/MEM MSc/1thesis/dems/CB18.tif")

# horizontal extent in metres
e <- ext(dem_check)
width_deg  <- e$xmax - e$xmin
height_deg <- e$ymax - e$ymin

# at -14.7 latitude: 1 degree longitude ≈ 107,000m, 1 degree latitude ≈ 111,000m
width_m  <- width_deg  * 107000
height_m <- height_deg * 111000

# pixel size in metres
pixel_size_m <- (width_m / ncol(dem_check) + height_m / nrow(dem_check)) / 2

# vertical range
relief <- global(dem_check, "max", na.rm=TRUE)$max - 
  global(dem_check, "min", na.rm=TRUE)$min

cat("Pixel size (m):", pixel_size_m, "\n")
cat("Total relief (m):", relief, "\n")
cat("Correct zscale:", pixel_size_m, "\n")




###plot per logger
dem <- prepare_dem("CB", 18)
hobo_sf <- load_hobo_locations("CB", 18, geom)
hobo_wgs84 <- st_transform(hobo_sf, "EPSG:4326")

png("logger_check.png", width=2400, height=2400)
par(mfrow=c(6,5), mar=c(1,1,2,1))

s <- ray_shade(dem$matrix, sunaltitude=45, sunangle=-70,
               zscale=0.004, maxsearch=200)
sr <- rast(s)
ext(sr) <- ext(dem$rast)
crs(sr) <- crs(dem$rast)

for(i in 1:nrow(hobo_wgs84)){
  plot(sr, legend=FALSE, axes=FALSE,
       main=hobo_sf$hoboID[i])
  plot(vect(hobo_wgs84[i,]), add=TRUE, 
       col="red", cex=3, pch=16)
}
dev.off()


# all 11 timesteps from your pipeline output
timesteps <- list(
  list(time="2018-12-04 11:00", alt=72.1,  az=-62.8),
  list(time="2018-12-04 12:00", alt=82.2,  az=-15.1),
  list(time="2018-12-04 13:00", alt=75.7,  az=56.4),
  list(time="2018-12-04 14:00", alt=62.6,  az=69.9),
  list(time="2018-12-04 15:00", alt=48.8,  az=72.9),
  list(time="2018-12-04 16:00", alt=35.0,  az=72.9),
  list(time="2018-12-04 17:00", alt=21.1,  az=71.5),
  list(time="2018-12-05 07:00", alt=17.1,  az=-70.8),
  list(time="2018-12-05 08:00", alt=30.8,  az=-72.5),
  list(time="2018-12-05 09:00", alt=44.7,  az=-73.0),
  list(time="2018-12-05 10:00", alt=58.5,  az=-71.1)
)

png("shadow_multipanel_labelled.png", width=3200, height=2400, res=150)
par(mfrow=c(3, 4), mar=c(1, 1, 3, 1))

for(ts in timesteps) {
  shadow <- ray_shade(dem$matrix,
                      sunaltitude = ts$alt,
                      sunangle    = ts$az,
                      zscale      = 0.005,
                      maxsearch   = 200)
  
  shadow_rast <- rast(shadow)
  ext(shadow_rast) <- ext(dem$rast)
  crs(shadow_rast) <- crs(dem$rast)
  
  plot(shadow_rast,
       main = paste0(ts$time, "\nalt=", ts$alt, "° az=", ts$az, "°"),
       legend = FALSE,
       axes = FALSE)
  
  plot(vect(hobo_wgs84), add=TRUE, col="red", cex=1.2, pch=16)
  
  text(st_coordinates(hobo_wgs84)[,1],
       st_coordinates(hobo_wgs84)[,2],
       labels = hobo_sf$hoboID,
       col = "white", cex = 0.55, font = 2, pos = 3)
}

# empty last panel for legend
plot.new()
legend("center",
       legend = c("Logger location"),
       col = "red", pch = 16, pt.cex = 1.5,
       bty = "n", cex = 1.2)

dev.off()
message("Saved shadow_multipanel_labelled.png")


