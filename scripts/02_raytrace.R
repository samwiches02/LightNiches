# ══════════════════════════════════════════════════════════════════════════════
#  02_raytrace.R — ray-trace every site once and build the cache table
#
#  Reads : output/hourly_inputs.rds, output/geom.rds, data/dems/<SITEyy>.tif
#  Writes: output/dems/<SITEyy>.rds      prepared DEM + ambient layer (for figures)
#          output/raytrace/<SITEyy>.rds  shading at every logger, every hour
#          output/cache.rds              hourly inputs + shading + measurement
#
#  This is the slow stage. Each site is saved as soon as it finishes, and a
#  site that is already saved is skipped (set FORCE_RAYTRACE <- TRUE in
#  config.R to redo). The measured light plays no part here.
#
#  [YOURS, EDITED] overall: this is your run_pipeline() outer and inner loops,
#  minus the Kd and prediction steps (now in 03). Changes are listed below.
#
#  PROVENANCE TAGS:  [YOURS] your code from master_script.R, logic unchanged
#                    [YOURS, EDITED] your code with listed changes   [NEW] Claude
#  Comments with no "[Claude]" marker are yours, copied from master_script.R.
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_raytrace.R"))
suppressPackageStartupMessages({ library(dplyr) })

stage_banner("02  RAY-TRACING")

# [NEW] Measure how ray_shade's sunangle maps onto compass bearings on our
#   grids, check it at four more angles, and stop if it doesn't hold. Every
#   shadow below is cast with this calibration. Saved for the record.
# [NEW] (1) how rayshader lays out its input, (2) a round-trip check, (3) how it
#   lays out its OUTPUT (found transposed on the CB test), (4) the sun direction.
sun_cal <- init_geometry_checks()
write.csv(data.frame(s = sun_cal$s, offset_deg = sun_cal$c),
          out_path("sun_calibration.csv"), row.names = FALSE)

hourly <- readRDS(out_path("hourly_inputs.rds"))
geom   <- readRDS(out_path("geom.rds"))
buffers <- c(BUFFER_M, BUFFER_SENS_M)

# ── raytrace_site ─────────────────────────────────────────────────────────────
# [YOURS, EDITED] — your run_pipeline() loops for one site/year
#   - [Claude] saves the prepared DEM and the result to disk (was dem_store
#     and results in memory)
#   - [Claude] sun at mid-hour; hours below SUN_ALT_MIN_DEG skipped here
#   - [Claude] also ray-traces the unrefracted sun (diagnostic) and extracts at
#     the sensitivity buffers
#   - [Claude] ambient extracted once per site instead of every hour (it does
#     not depend on the sun)
#
#outer loop: load spatial data once per site/year
#when they are outside hourly loop we only do this once, since it takes a lot of time (it can hang R if done inside the loop)
raytrace_site <- function(key) {
  t0 <- Sys.time()
  dem <- prepare_dem(key)
  save_dem(dem, out_path("dems", paste0(key, ".rds")))

  hobo     <- load_hobo_locations(geom, key)
  # [NEW] positions and DEM must be in the same place: the reef disc is centred
  #   on the DEM, so count loggers outside it. A whole site outside means the
  #   positions and the DEM are in different coordinate frames (likely for a
  #   batch delivered in another projection).
  e   <- terra::ext(dem$rast)
  xy  <- sf::st_coordinates(hobo)
  off <- sqrt((xy[, 1] - (e$xmin + e$xmax) / 2)^2 + (xy[, 2] - (e$ymin + e$ymax) / 2)^2)
  n_off <- sum(off > REEF_RADIUS_M)
  message(sprintf("  logger distance from reef-disc centre: median %.1f m, max %.1f m (disc radius %.1f m)",
                  median(off), max(off), REEF_RADIUS_M))
  if (n_off == length(off))
    stop(key, ": none of the ", n_off, " logger positions fall on the reef disc. The ",
         "positions in geom_new.csv and the DEM are in different coordinate frames.",
         call. = FALSE)
  if (n_off > 0) warning(key, ": ", n_off, " of ", length(off),
                         " logger positions fall outside the reef disc", call. = FALSE)
  centroid <- get_site_centroid(dem$rast)
  amb      <- dem$ambient; amb[!dem$reef_mask] <- NA

  # ambient is sun-independent: extract once per buffer
  amb_cols <- lapply(buffers, function(b) extract_at_loggers(amb, dem$rast, hobo, b))

  #get list of hours to raytrace, we only need one shadow map per hour
  hours <- sort(unique(hourly$datetime_pos[hourly$sitekey == key]))
  rows  <- list()
  #inner loop, one raytracing per hour
  #raytraces once and extracts shadow + ambient values for all loggers from
  #that single shadow map (one raytrace, many extractions).
  for (i in seq_along(hours)) {
    sun <- sun_at_midhour(hours[i], centroid$lat, centroid$lon)
    alt_deg <- sun$altitude * 180 / pi
    if (alt_deg <= SUN_ALT_MIN_DEG) next

    s_refr <- make_shadow(dem, sun, refract = TRUE)
    s_air  <- if (RUN_UNREFRACTED_DIAGNOSTIC) make_shadow(dem, sun, refract = FALSE)

    df <- data.frame(sitekey = key, hoboID_env = hobo$hoboID,
                     datetime_pos = hours[i],
                     sun_altitude = sun$altitude, sun_azimuth = sun$azimuth)
    for (j in seq_along(buffers)) {
      sfx <- if (j == 1) "" else buffer_suffix(buffers[j])
      df[[paste0("shadow_refr", sfx)]] <- extract_at_loggers(s_refr, dem$rast, hobo, buffers[j])
      df[[paste0("ambient",     sfx)]] <- amb_cols[[j]]
      if (RUN_UNREFRACTED_DIAGNOSTIC && j == 1)
        df$shadow_air <- extract_at_loggers(s_air, dem$rast, hobo, buffers[j])
    }
    rows[[length(rows) + 1]] <- df
    message(sprintf("  %s %s  alt %.1f° (in water %.1f°)", key,
                    format(hours[i], "%m-%d %H:%M"), alt_deg,
                    raytrace_angles(sun)$alt_deg))
  }
  out <- bind_rows(rows)
  if (!RUN_UNREFRACTED_DIAGNOSTIC) out$shadow_air <- NA_real_
  saveRDS(out, out_path("raytrace", paste0(key, ".rds")))
  message(sprintf("  %s done: %d hours, %.1f min", key, length(rows),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  out
}

# [NEW] skip sites already saved, so an interrupted run resumes
for (key in SITES) {
  f <- out_path("raytrace", paste0(key, ".rds"))
  if (file.exists(f) && !FORCE_RAYTRACE) { message("\n", key, ": cached, skipping"); next }
  message("\n── ", key, " ──")
  raytrace_site(key)
}

# ── Build the cache: inputs + shading + measured light, one row per logger-hour
# [YOURS, EDITED] — your VALIDATION JOIN, keyed on datetime + logger + sitekey
#   so loggers can't collide across sites. [Claude] Reference loggers are NOT
#   dropped here; 03 decides which rows are evaluated.
shade <- bind_rows(lapply(SITES, function(k) readRDS(out_path("raytrace", paste0(k, ".rds")))))
cache <- hourly %>%
  inner_join(shade, by = c("sitekey", "hoboID_env", "datetime_pos"))

dropped <- hourly %>% distinct(sitekey, hoboID_env) %>%
  anti_join(distinct(cache, sitekey, hoboID_env), by = c("sitekey", "hoboID_env"))
message("\nCache: ", nrow(cache), " logger-hours, ",
        nrow(distinct(cache, sitekey, hoboID_env)), " loggers. Loggers dropped (no position): ",
        nrow(dropped))
saveRDS(cache, out_path("cache.rds"))

write_manifest("02_raytrace", extra = c(
  sprintf("cache rows: %d", nrow(cache)),
  sprintf("loggers dropped for missing positions: %d", nrow(dropped)),
  sprintf("sun calibration: s = %+d, offset = %.0f deg (checked at 45/135/225/315)",
          sun_cal$s, sun_cal$c)))
message("02 done.")
