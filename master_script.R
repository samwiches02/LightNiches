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

KD_CONSTANT <- 0.25 #KD Constant -universal diffuse attenutation coefficient (m⁻¹) (still need to look at jervlov water numbers)
Z_SCALE <- 0.02 #the zscale
TIMESTEP_INTERVAL <- "1 hour" #the timestep interval to measure
DOWNSAMPLE_FACTOR <- 10 #by how many times to downsample the raster for faster processing speed, up to 4 does not severely affect quality, 3 makes it 9x faster

USE_AMBIENT <- FALSE #toggle for ambient occlusion layer, if true, slower processing speed
MAXSEARCH_FACT <- 100 #change maxsearch distance, default 200, basically, every pixel shoots up towards the sun at x distance, and if it hits something within that distance, it is a shadow, 
SHADOW_GAMMA <- 1 #shadow gamma which acts to darken areas by applying a power of gamma to the shadow*ambient layer,this will crush the shadow layer (hopefully) to make darker stronger and does not affect lighter areas as much




# ══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# reusable code chunk recipes that I made to be called later 
# ══════════════════════════════════════════════════════════════════════════════
# ── ZSCALE CALCULATOR ─────────────────────────────────────────────────────────
calc_zscale <- function(dem_path, downsample_factor) {
  dem        <- rast(dem_path)
  e          <- ext(dem)
  width_m    <- (e$xmax - e$xmin) * 107000
  height_m   <- (e$ymax - e$ymin) * 111000
  pixel_size <- (width_m / ncol(dem) + height_m / nrow(dem)) / 2
  zscale     <- pixel_size * downsample_factor
  
  cat("Native pixel size (m):       ", round(pixel_size, 6), "\n")
  cat("Downsample factor:           ", downsample_factor, "\n")
  cat("Downsampled pixel size (m):  ", round(pixel_size * downsample_factor, 6), "\n")
  cat("Recommended Z_SCALE:         ", round(zscale, 6), "\n")
  
  invisible(zscale)
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
  
  # ── NEW: circular trim then NA fill ──────────────────────────────────────
  # Step 1: apply circular mask to remove border artefacts
  trim_to_circle <- function(mat, shrink = 0.90) {
    nrows      <- nrow(mat)
    ncols      <- ncol(mat)
    center_row <- nrows / 2
    center_col <- ncols / 2
    radius     <- min(center_row, center_col) * shrink  # shrink inward by 10%
    row_idx    <- matrix(rep(1:nrows, ncols), nrow = nrows, ncol = ncols)
    col_idx    <- matrix(rep(1:ncols, nrows), nrow = nrows, ncol = ncols, byrow = TRUE)
    distances  <- sqrt((row_idx - center_row)^2 + (col_idx - center_col)^2)
    mat[distances > radius] <- NA
    mat
  }
  
  dem_matrix <- trim_to_circle(dem_matrix, shrink = 0.90)
  
  # Step 2: fill remaining NAs with median so rayshader sees flat continuation
  # now that border artefacts are removed, median is a safe neutral value
  dem_matrix[is.na(dem_matrix)] <- median(dem_matrix, na.rm = TRUE)
  
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
make_shadow <- function(dem_matrix, sun, zscale = Z_SCALE) {
  if (sun$altitude <= 0.1) {
    message("    (sun below horizon — skipping raytrace)")
    return(NULL)
  }
  
  shadow <- ray_shade(
    dem_matrix,
    sunaltitude = sun$altitude * 180 / pi, #gets sun altitude in degrees from radians
    sunangle    = sun$azimuth  * 180 / pi, #gets sun azimuth in degrees from radians
    zscale      = zscale,          
    maxsearch   = MAXSEARCH_FACT #limits how each ray travels, so smaller features have more importance (but this takes longer)
  )
  
  if (USE_AMBIENT) {
    ambient <- ambient_shade(dem_matrix, zscale = zscale)  #calculates a new shadow layer based on ambient shade, making more enclosed spots darker regardless of sun angl
    shadow  <- (shadow * ambient) ^ SHADOW_GAMMA  #darkens areas by applying a power of gamma to the shadow*ambient layer,this will crush the shadow layer (hopefully) to make darker stronger and lighter lighter
  }
  
  message("    shadow range: ", round(min(shadow, na.rm = TRUE), 3),
          " to ",               round(max(shadow, na.rm = TRUE), 3),
          " | mean: ",          round(mean(shadow, na.rm = TRUE), 3))
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
  crs(shadow_rast) <- crs(dem_rast)
  
  hobo_wgs84 <- st_transform(hobo_sf, "EPSG:4326")
  
  # Buffer scaled to match physical sensor casing size (~5cm radius)
  # regardless of what downsample factor was used.
  # At -14.7° latitude, 1 degree ≈ 111,000m
  # So 0.05m (5cm) = 0.05 / 111000 degrees
  buffer_deg <- 0.05 / 111000
  
  hobo_buffered <- st_buffer(hobo_wgs84, dist = buffer_deg)
  
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

#HOURLY THING TO CUT DOWN EVEN MORE NOISE 
#env_long_test <- env_long %>%
#  filter(site == "HS", year == 17) %>%
#  filter(datetime_pos >= as.POSIXct("2017-11-17 11:00:00", tz = "Australia/Brisbane"),
#         datetime_pos <  as.POSIXct("2017-11-17 18:00:00", tz = "Australia/Brisbane"))

# ── FILTER DATA ───────────────────────────────────────────────────────────────
env_long_test <- env_long %>%
  filter(site == "TM", year == 19)

# derived immediately from filter — never hardcode these
active_site <- unique(env_long_test$site)
active_year <- unique(env_long_test$year)

# ── DEM VERIFICATION ──────────────────────────────────────────────────────────
calc_zscale(
  dem_path          = file.path(DEM_DIR, paste0(active_site, active_year, ".tif")),
  downsample_factor = DOWNSAMPLE_FACTOR
)

dem_check   <- prepare_dem(active_site, active_year)
real_pixels <- dem_check$matrix[dem_check$matrix != median(dem_check$matrix, na.rm = TRUE)]

hist(real_pixels, breaks = 100,
     main   = paste0(active_site, active_year, " DEM elevation distribution (real pixels only)"),
     xlab   = "Elevation (m)", col = "grey80", border = "white")

cat("\nElevation quantiles:\n")
print(quantile(real_pixels, probs = c(0.01, 0.25, 0.50, 0.75, 0.99), na.rm = TRUE))
cat("Total relief:", round(diff(range(real_pixels, na.rm = TRUE)), 3), "m\n")

# ── RUN PIPELINE ──────────────────────────────────────────────────────────────
output <- run_pipeline(env_long_test, geom)

# ── VALIDATION JOIN ───────────────────────────────────────────────────────────
validation <- output %>%
  left_join(
    env_long_test %>%
      mutate(hoboID_env = as.character(hoboID)) %>%
      mutate(datetime_pos = floor_date(datetime_pos, TIMESTEP_INTERVAL)) %>%
      group_by(datetime_pos, hoboID_env) %>%
      summarise(light_mol = mean(light_mol, na.rm = TRUE), .groups = "drop"),
    by = c("datetime_pos", "hoboID_env")
  ) %>%
  filter(!is.na(pred_par), !is.na(light_mol))

#check if rows filled
nrow(validation)




# ── PLOTS ─────────────────────────────────────────────────────────────────────
# 1. time series per logger
# 2. shadow multipanel
# 3. DLI scatter
# 4. time-of-day scatter (fix hardcoded title here!!!
# ══════════════════════════════════════════════════════════════════════════════
#  PLOT — measured vs predicted PAR per logger
# ══════════════════════════════════════════════════════════════════════════════
dev.off()
dev.new()

site_name  <- unique(validation$site)
date_range <- format(range(validation$datetime_pos), "%Y-%m-%d %H:%M")
plot_title <- paste0(
  site_name, " | ", date_range[1], " to ", date_range[2], "\n",
  "Kd=", KD_CONSTANT,
  " | gamma=", SHADOW_GAMMA,
  " | zscale=", Z_SCALE,
  " | interval=", TIMESTEP_INTERVAL,
  " | ambient layer=", USE_AMBIENT,
  " | downsampled by x=", DOWNSAMPLE_FACTOR
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



# ══════════════════════════════════════════════════════════════════════════════
#  SHADOWMAPS + CHART OF EACH POINT IN TIME FROM PIPELINE (uses correct sun angles)
# ══════════════════════════════════════════════════════════════════════════════
sun_steps <- output %>%
  distinct(datetime_pos, sun_altitude, sun_azimuth) %>%
  arrange(datetime_pos)

n_steps <- nrow(sun_steps)
n_cols  <- 4
n_rows  <- ceiling(n_steps / n_cols)

draw_shadow_multipanel <- function() {
  layout_mat <- matrix(0, nrow = n_rows * 2 + 2, ncol = n_cols)
  layout_mat[1, ] <- 1
  panel <- 2
  for (r in seq_len(n_rows)) {
    for (c in seq_len(n_cols)) {
      step_i <- (r - 1) * n_cols + c
      if (step_i <= n_steps) {
        layout_mat[r * 2,     c] <- panel
        layout_mat[r * 2 + 1, c] <- panel + 1
        panel <- panel + 2
      }
    }
  }
  legend_cell <- panel
  layout_mat[n_rows * 2 + 2, ] <- legend_cell
  layout(layout_mat, heights = c(0.6, rep(c(5, 1.5), n_rows), 1.5))
  
  # ── title panel ──
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5,
       paste0(unique(output$site[!is.na(output$site)]),
              unique(output$year[!is.na(output$year)]),
              "  |  Kd=", KD_CONSTANT,
              "  |  zscale=", Z_SCALE,
              "  |  downsample=", DOWNSAMPLE_FACTOR,
              "  |  maxsearch=", MAXSEARCH_FACT,
              "  |  ambient=", USE_AMBIENT,
              "  |  gamma=", SHADOW_GAMMA),
       cex = 1.2, font = 2, adj = 0.5)
  
  for (i in seq_len(n_steps)) {
    ts      <- sun_steps[i, ]
    alt_deg <- ts$sun_altitude * 180 / pi
    az_deg  <- ts$sun_azimuth  * 180 / pi
    
    # ── shadow map panel ──
    par(mar = c(0, 1, 5, 1))
    shadow <- ray_shade(dem$matrix,
                        sunaltitude = alt_deg,
                        sunangle    = az_deg,
                        zscale      = Z_SCALE,
                        maxsearch   = MAXSEARCH_FACT)
    
    shadow_rast <- rast(shadow)
    ext(shadow_rast) <- ext(dem$rast)
    crs(shadow_rast) <- crs(dem$rast)
    
    plot(shadow_rast,
         main     = paste0(format(ts$datetime_pos, "%Y-%m-%d %H:%M"),
                           "\nalt=", round(alt_deg, 1), "° az=", round(az_deg, 1), "°"),
         legend   = FALSE, axes = FALSE, cex.main = 1.0)
    
    plot(vect(hobo_wgs84), add = TRUE,
         col = "red", bg = NA, cex = 1.5, pch = 21, lwd = 2)
    
    text(st_coordinates(hobo_wgs84)[, 1],
         st_coordinates(hobo_wgs84)[, 2],
         labels = hobo_sf$hoboID,
         col = "yellow", cex = 0.9, font = 2, pos = 3)
    
    # ── bar chart panel ──
    par(mar = c(4, 3, 0.5, 1))
    
    step_data <- validation %>%
      filter(datetime_pos == ts$datetime_pos) %>%
      arrange(hoboID_env)
    
    if (nrow(step_data) > 0) {
      bar_mat <- rbind(step_data$pred_par, step_data$light_mol)
      bar_mat[is.na(bar_mat)] <- 0
      
      barplot(bar_mat,
              beside    = TRUE,
              col       = c("red", "blue"),
              names.arg = step_data$hoboID_env,
              cex.names = 0.8,
              cex.axis  = 0.9,
              las       = 2,
              ylim      = c(0, max(c(validation$pred_par,
                                     validation$light_mol), na.rm = TRUE)),
              border    = NA)
    } else {
      plot.new()
    }
  }
  
  # ── legend panel ──
  par(mar = c(2, 4, 2, 4))
  plot.new()
  
  n_grad  <- 100
  grad    <- hcl.colors(n_grad, palette = "viridis")
  x_left  <- 0.05
  x_right <- 0.55
  y_bot   <- 0.3
  y_top   <- 0.65
  x_steps <- seq(x_left, x_right, length.out = n_grad + 1)
  
  for (j in seq_len(n_grad)) {
    rect(x_steps[j], y_bot, x_steps[j + 1], y_top,
         col = grad[j], border = NA)
  }
  rect(x_left, y_bot, x_right, y_top, border = "black", lwd = 0.5)
  text(x_left,                  y_bot - 0.15, "0\n(full shadow)", cex = 1.1, adj = 0.5)
  text(x_right,                 y_bot - 0.15, "1\n(full sun)",    cex = 1.1, adj = 0.5)
  text((x_left + x_right) / 2, y_top + 0.15, "Shadow fraction",  cex = 1.3, font = 2)
  
  legend(0.65, 0.8,
         legend = c("Predicted PAR", "Measured PAR"),
         fill   = c("red", "blue"),
         border = NA,
         bty    = "n",
         cex    = 1.2,
         title  = "Bar chart key")
  
  par(mfrow = c(1, 1))
}

# draw to screen
dev.off()
dev.new()
draw_shadow_multipanel()

# save to file
png("shadow_multipanel_pipeline.png",
    width  = n_cols * 800,
    height = n_rows * 1100 + 200,
    res    = 150)
draw_shadow_multipanel()
dev.off()
message("Saved shadow_multipanel_pipeline.png")



# ── DLI PLOT ─────────────────────────────────────────────────────────────────
# pull site and year automatically from whatever is in env_long_test
active_site <- unique(env_long_test$site)
active_year <- unique(env_long_test$year)

dev.off()
dev.new()
par(mar = c(5, 5, 4, 2))

plot_title <- paste0(
  "Daily light integral, measured vs predicted, ",
  paste(active_site, active_year, sep = "", collapse = ", "), "\n",
  "Kd=", KD_CONSTANT,
  " | zscale=", Z_SCALE,
  " | downsample=", DOWNSAMPLE_FACTOR,
  " | gamma=", SHADOW_GAMMA,
  " | ambient=", USE_AMBIENT,
  " | maxsearch=", MAXSEARCH_FACT
)

# measured DLI: sum of 5-min µmol m-2 s-1 readings × 300s per interval ÷ 1,000,000
dli_measured <- aggregate(
  light_mol ~ hoboID + as.Date(datetime_pos),
  data = env_long[env_long$site %in% active_site & 
                    env_long$year %in% active_year, ],
  FUN = function(x) sum(x, na.rm = TRUE) * 300 / 1e6
)
colnames(dli_measured) <- c("hoboID", "date", "dli_measured")

# predicted DLI: sum of hourly µmol m-2 s-1 readings × 3600s per interval ÷ 1,000,000
validation$date <- as.Date(validation$datetime_pos)
dli_predicted <- aggregate(
  pred_par ~ hoboID_env + date,
  data = validation,
  FUN = function(x) sum(x, na.rm = TRUE) * 3600 / 1e6
)
colnames(dli_predicted) <- c("hoboID", "date", "dli_predicted")

# join and plot
dli <- merge(dli_measured, dli_predicted, by = c("hoboID", "date"))

plot(dli$dli_measured, dli$dli_predicted,
     pch  = 16, col  = "steelblue",
     xlab = "Measured DLI (mol m-2 d-1)",
     ylab = "Predicted DLI (mol m-2 d-1)",
     main = plot_title, cex.main = 0.85)

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
     main = paste0("Predicted vs Observed PAR by time of day, ",
                   active_site, active_year,
                   " | Kd=", KD_CONSTANT,
                   " | zscale=", Z_SCALE))

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



