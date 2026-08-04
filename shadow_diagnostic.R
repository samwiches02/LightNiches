# ══════════════════════════════════════════════════════════════════════════════
#  SINGLE-SENSOR DIAGNOSTIC
#  Pick one random logger, trace every component of pred_par through the day,
#  and zoom into the exact pixels it samples from the shadow map.
#
#  RUN THIS AFTER the pipeline (needs: output, validation, dem_store,
#  env_long_test, and the pipeline constants BUFFER_M / W_DIRECT / W_AMBIENT /
#  USE_AMBIENT / MAXSEARCH_PX). Nothing here is written back into the pipeline.
# ══════════════════════════════════════════════════════════════════════════════


# ── pick a specific sensor ───────────────────────────────────────────────────
pick_key <- "TB19"      # sitekey: "CB19","GT19","HS19","N319","TB19","TM19"
pick_id  <- "H50"       # logger ID at that site

# ── manual position override ──────────────────────────────────────────────────
# If section 6 told you sensor X wants to shift +0.15 m east, -0.10 m north,
# set MANUAL_OFFSET_M = c(0.15, -0.10) here to see what the trace/maps/sweep
# look like AS IF the sensor sat there. Set to c(0, 0) to keep the original.
MANUAL_OFFSET_M <- c(dx = 0.00, dy = 0.00)   # metres in UTM (+east, +north)

if (any(MANUAL_OFFSET_M != 0)) {
  cat(sprintf("\n** Applying manual offset to %s: dx=%+.3f m, dy=%+.3f m **\n",
              pick_id, MANUAL_OFFSET_M["dx"], MANUAL_OFFSET_M["dy"]))
  
  # shift the picked sensor in the UTM-projected hobo layer
  idx <- which(hobo_utm$hoboID == pick_id)
  orig_geom <- st_geometry(hobo_utm)[[idx]]
  new_coords <- st_coordinates(orig_geom) + matrix(MANUAL_OFFSET_M, nrow = 1)
  st_geometry(hobo_utm)[[idx]] <- st_point(as.numeric(new_coords))
  
  # also shift the WGS84 copy used by some plotting blocks (project the new pt)
  shifted_wgs <- st_transform(st_sfc(st_point(as.numeric(new_coords)),
                                     crs = UTM_CRS), 4326)
  st_geometry(hobo_sf)[[idx]] <- st_geometry(shifted_wgs)[[1]]
}

stopifnot(pick_key %in% validation$sitekey)
stopifnot(pick_id %in% validation$hoboID_env[validation$sitekey == pick_key])

# to force a specific logger instead, hardcode e.g.:
# pick_key <- "TB18"; pick_id <- "H34"

# ── pull this site's cached spatial objects from dem_store ────────────────────
store     <- dem_store[[pick_key]]
dem       <- store$dem          # $matrix, $rast, $zscale, $reef_mask, $ambient
hobo_sf   <- store$hobo_sf
hobo_utm  <- store$hobo_utm     # already in UTM 55S

# ── read the pipeline's own constants so the diagnostic can never drift ───────
# (these are globals defined in master_script.R; we just reference them)
BUFFER_R  <- BUFFER_M           # pipeline buffer radius (master uses st_buffer dist = BUFFER_M)


# pipeline uses MAXSEARCH_M (metres) and converts to pixels per-DEM via zscale.
# Derive the pixel-equivalent for THIS site so every ray_shade call below
# uses the same physical ray distance the pipeline used.
MAXSEARCH_PX <- max(1, round(MAXSEARCH_M / dem$zscale))
cat("MAXSEARCH_M =", MAXSEARCH_M, "m -> ", MAXSEARCH_PX, "px on this DEM\n")
# the pipeline's effective-light rule, replicated EXACTLY (honours USE_AMBIENT).
# shadow + ambient are vectors of the same length; returns blended eff_light.
eff_light_fn <- function(shadow_frac, ambient_frac) {
  if (USE_AMBIENT) {
    W_DIRECT * shadow_frac + W_AMBIENT * ambient_frac
  } else {
    shadow_frac                 # ambient off: direct beam only
  }
}

cat("\n=====================================================\n")
cat("DIAGNOSTIC SENSOR:", pick_id, " @ site/year:", pick_key, "\n")
cat("pipeline buffer radius:", BUFFER_R, "m  |  USE_AMBIENT =", USE_AMBIENT,
    " (w_dir=", W_DIRECT, " w_amb=", W_AMBIENT, ")\n", sep = "")
cat("=====================================================\n")

# ── 1. component table for the full day cycle ─────────────────────────────────
# validation already holds everything the prediction is built from:
#   surface_par      – top-of-water PAR for that hour
#   shadow_fraction  – direct-beam (ray_shade) value sampled under the sensor
#   ambient_fraction – diffuse (ambient_shade) value sampled under the sensor
#   eff_light        – W_DIRECT*shadow + W_AMBIENT*ambient (or shadow only)
#   water_atten      – exp(-Kd * depth)
#   pred_par         = surface_par * eff_light * water_atten
#   light_mol        – the MEASURED value to compare against
trace_tbl <- validation %>%
  filter(sitekey == pick_key, hoboID_env == pick_id) %>%
  arrange(datetime_pos) %>%
  transmute(
    time          = format(datetime_pos, "%m-%d %H:%M"),
    sun_alt_deg   = round(sun_altitude * 180 / pi, 1),
    surface_par   = round(surface_par, 1),
    Kd            = round(Kd, 4),
    depth_m       = round(water_depth_m, 2),
    water_atten   = round(water_atten, 4),        # exp(-Kd*depth)
    shadow_frac   = round(shadow_fraction, 4),    # direct-beam pixel value(s)
    ambient_frac  = round(ambient_fraction, 4),   # diffuse pixel value(s)
    eff_light     = round(eff_light_fn(shadow_fraction, ambient_fraction), 4),
    pred_par      = round(pred_par, 1),
    measured      = round(light_mol, 1),
    ratio_p_over_m = round(pred_par / light_mol, 2)
  )

cat("\nComponent trace (pred = surface_par x eff_light x water_atten),\n")
if (USE_AMBIENT) {
  cat("  eff_light = ", W_DIRECT, "*shadow + ", W_AMBIENT, "*ambient :\n\n", sep = "")
} else {
  cat("  eff_light = shadow only (USE_AMBIENT = FALSE) :\n\n", sep = "")
}
print(as.data.frame(trace_tbl), row.names = FALSE)

cat("\n-- day summary --\n")
cat("mean shadow_frac :", round(mean(trace_tbl$shadow_frac,  na.rm=TRUE), 3), "\n")
cat("mean ambient_frac:", round(mean(trace_tbl$ambient_frac, na.rm=TRUE), 3), "\n")
cat("mean eff_light   :", round(mean(trace_tbl$eff_light,    na.rm=TRUE), 3), "\n")
cat("mean water_atten :", round(mean(trace_tbl$water_atten,  na.rm=TRUE), 3), "\n")
cat("depth (m)        :", round(mean(trace_tbl$depth_m,      na.rm=TRUE), 2), "\n")
cat("pred DLI  (mol)  :", round(sum(trace_tbl$pred_par,  na.rm=TRUE) * 3600 / 1e6, 2), "\n")
cat("meas DLI  (mol)  :", round(sum(trace_tbl$measured,  na.rm=TRUE) * 3600 / 1e6, 2), "\n")


# ── 2. pixel-level zoom: what exactly does the sensor sample? ──────────────────
# Recreate the shadow raster for ONE chosen hour, then pull the pixel window
# around the sensor so you can see the raw values that get averaged in the buffer.

# choose the timestep to inspect — default = solar noon (highest sun)
# output is pooled over ALL site/years, so filter to THIS sensor's site first.
sun_steps <- output %>%
  filter(sitekey == pick_key) %>%
  distinct(datetime_pos, sun_altitude, sun_azimuth) %>%
  arrange(datetime_pos)
inspect_row <- sun_steps[which.max(sun_steps$sun_altitude), ]   # or pick by index
inspect_time <- inspect_row$datetime_pos
alt_deg <- inspect_row$sun_altitude * 180 / pi
az_deg  <- inspect_row$sun_azimuth  * 180 / pi

cat("\n=====================================================\n")
cat("PIXEL ZOOM at", format(inspect_time, "%Y-%m-%d %H:%M"),
    " alt=", round(alt_deg,1), " az=", round(az_deg,1), "\n")
cat("=====================================================\n")

# rebuild the shadow matrix for this hour (same call the pipeline uses)
shadow <- ray_shade(dem$matrix,
                    sunaltitude = alt_deg,
                    sunangle    = az_deg,
                    zscale      = dem$zscale,
                    maxsearch   = MAXSEARCH_PX)

# wrap into a georeferenced raster EXACTLY as extract_hobo_shadow does
# (note: your current pipeline uses rast(t(shadow)) — match that here)
shadow_rast <- rast(t(shadow))
ext(shadow_rast) <- ext(dem$rast)
crs(shadow_rast) <- crs(dem$rast)            # UTM 55S

# sensor location in UTM (same CRS as the shadow raster)
this_hobo <- hobo_utm[hobo_utm$hoboID == pick_id, ]
xy <- st_coordinates(this_hobo)
cat("Sensor UTM x/y:", round(xy[1], 2), round(xy[2], 2), "\n")

# the cell the sensor centre falls in
cell <- terra::cellFromXY(shadow_rast, xy)
centre_val <- shadow_rast[cell][1, 1]
cat("Centre pixel shadow value:", round(centre_val, 4), "\n")

# pixel resolution (metres) so you can see how big the buffer is vs a pixel
px <- terra::res(shadow_rast)
cat("Pixel size (m):", round(px[1], 3), "x", round(px[2], 3),
    "  | buffer radius =", BUFFER_R, "m\n")
cat("  (if buffer << pixel, the buffer-mean just equals the single centre pixel)\n")

# the buffered-mean value the pipeline actually assigns (BUFFER_R = pipeline buffer):
hobo_buf <- st_buffer(this_hobo, dist = BUFFER_R)
buf_val  <- terra::extract(shadow_rast, vect(hobo_buf), fun = mean, na.rm = TRUE)[,2]
cat("Buffer-mean shadow value (what pipeline uses):", round(buf_val, 4), "\n")

# ── 2b. print a small window of raw pixels around the sensor ──────────────────
rc <- terra::rowColFromCell(shadow_rast, cell)
r0 <- rc[1]; c0 <- rc[2]
half <- 4                                          # 9x9 window
rows <- (r0 - half):(r0 + half)
cols <- (c0 - half):(c0 + half)
rows <- rows[rows >= 1 & rows <= nrow(shadow_rast)]
cols <- cols[cols >= 1 & cols <= ncol(shadow_rast)]

window <- matrix(NA_real_, length(rows), length(cols),
                 dimnames = list(paste0("r", rows), paste0("c", cols)))
for (ii in seq_along(rows))
  for (jj in seq_along(cols))
    window[ii, jj] <- shadow_rast[terra::cellFromRowCol(shadow_rast, rows[ii], cols[jj])][1,1]

cat("\nRaw shadow pixels around sensor (centre = r", r0, "c", c0, "):\n", sep="")
print(round(window, 3))


# ── 3. hourly measured vs predicted, fully labelled ───────────────────────────
dev.new(width = 8, height = 5.5)
par(mar = c(7, 4.5, 4, 2))
with(trace_tbl, {
  yr <- range(c(measured, pred_par), na.rm = TRUE)
  plot(seq_along(time), measured, type = "b", col = "blue", lwd = 2, pch = 16,
       ylim = yr, xaxt = "n", xlab = "", ylab = expression(PAR~(mu*mol~m^-2~s^-1)),
       main = paste0(pick_key, " / ", pick_id,
                     " — hourly measured vs predicted PAR\n",
                     "Kd=", round(Kd[1], 3), "  depth=", round(depth_m[1], 2), " m  |  ",
                     if (USE_AMBIENT)
                       paste0("eff = ", W_DIRECT, "*shadow + ", W_AMBIENT, "*ambient")
                     else "eff = shadow only",
                     "  |  buffer ", BUFFER_R, " m"))
  axis(1, at = seq_along(time), labels = time, las = 2, cex.axis = 0.65)
  lines(seq_along(time), pred_par, type = "b", col = "red", lwd = 2, lty = 2, pch = 1)
  
  # label each point with its value so nothing is unlabelled
  text(seq_along(time), measured, labels = round(measured), pos = 3,
       cex = 0.55, col = "blue")
  text(seq_along(time), pred_par, labels = round(pred_par), pos = 1,
       cex = 0.55, col = "red")
  
  # vertical guides showing the predicted-vs-measured gap each hour
  segments(seq_along(time), measured, seq_along(time), pred_par,
           col = "grey70", lty = 3)
  
  legend("topright", c("Measured", "Predicted", "hourly gap"),
         col = c("blue", "red", "grey70"), lty = c(1, 2, 3),
         pch = c(16, 1, NA), lwd = 2, bty = "n", cex = 0.8)
})


# ── 4. ROW of zoomed shadow maps, one per timestep ────────────────────────────
# rebuilds the shadow raster at every hour, crops a tight window around the
# sensor's buffer, and marks the sensor as an open circle.
# ZOOM_HALF sets the half-width of the window in metres. It is floored so the
# full buffer ring is always visible with a little margin around it.
buffer_r  <- BUFFER_R                             # matches pipeline buffer radius
ZOOM_HALF <- max(0.08, buffer_r * 1.6)            # never crops the ring off

# per-timestep measured + predicted for THIS sensor, keyed by datetime,
# so each panel title can show what the sensor actually recorded that hour
meas_lookup <- validation %>%
  filter(sitekey == pick_key, hoboID_env == pick_id) %>%
  select(datetime_pos, light_mol, pred_par)

n_steps <- nrow(sun_steps)
n_cols  <- min(n_steps, 4)
n_rows  <- ceiling(n_steps / n_cols)

dev.new(width = n_cols * 4.0, height = n_rows * 4.4 + 0.5)
par(mfrow = c(n_rows, n_cols), mar = c(1, 1, 5, 1), oma = c(0,0,2,0))

# circle helper for drawing the buffer footprint
circle_xy <- function(cx, cy, r, n = 60) {
  th <- seq(0, 2*pi, length.out = n)
  list(x = cx + r*cos(th), y = cy + r*sin(th))
}

for (k in seq_len(n_steps)) {
  trow    <- sun_steps[k, ]
  k_alt   <- trow$sun_altitude * 180 / pi
  k_az    <- trow$sun_azimuth  * 180 / pi
  
  k_shadow <- ray_shade(dem$matrix,
                        sunaltitude = k_alt,
                        sunangle    = k_az,
                        zscale      = dem$zscale,
                        maxsearch   = MAXSEARCH_PX)
  k_rast <- rast(t(k_shadow))
  ext(k_rast) <- ext(dem$rast)
  crs(k_rast) <- crs(dem$rast)
  
  # value the pipeline would assign this hour (buffer mean)
  k_val <- terra::extract(k_rast, vect(hobo_buf), fun = mean, na.rm = TRUE)[,2]
  
  # pull this sensor's measured + predicted PAR for this hour (may be absent)
  mp <- meas_lookup[meas_lookup$datetime_pos == trow$datetime_pos, ]
  meas_txt <- if (nrow(mp)) round(mp$light_mol[1]) else NA
  pred_txt <- if (nrow(mp)) round(mp$pred_par[1])  else NA
  
  win  <- ext(c(xy[1]-ZOOM_HALF, xy[1]+ZOOM_HALF,
                xy[2]-ZOOM_HALF, xy[2]+ZOOM_HALF))
  zoom <- crop(k_rast, win)
  
  plot(zoom, col = hcl.colors(100, "viridis"),
       range = c(0, 1), legend = FALSE, axes = FALSE,
       main = paste0(pick_key, " / ", pick_id, "\n",
                     format(trow$datetime_pos, "%m-%d %H:%M"),
                     "  alt=", round(k_alt, 1), "\u00b0\n",
                     "shadow=", round(k_val, 3),
                     "  meas=", meas_txt, "  pred=", pred_txt),
       cex.main = 0.95)
  
  # buffer footprint + sensor as an OPEN circle
  bc <- circle_xy(xy[1], xy[2], buffer_r)
  lines(bc$x, bc$y, col = "white", lwd = 1, lty = 3)        # the pipeline buffer
  points(xy[1], xy[2], pch = 1, col = "red", cex = 2, lwd = 2)  # open circle
}

mtext(paste0(pick_key, " / ", pick_id, " — shadow buffer zoom (+/-", ZOOM_HALF,
             " m), dotted = ", buffer_r*100, " cm pipeline buffer"),
      outer = TRUE, font = 2, cex = 0.9)

par(mfrow = c(1,1))

cat("\nDone. Edit `pick_key`/`pick_id` to lock a sensor, `ZOOM_HALF` to change\n")
cat("zoom level, or `inspect_row` to change the single-hour pixel-window printout.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  5. BUFFER-RADIUS SENSITIVITY SWEEP (for the chosen sensor)
#  How much does the shadow value — and therefore pred_par — depend on the
#  buffer size? Sweeps from a single centre pixel up to the pipeline buffer,
#  across the day. The pipeline buffer (BUFFER_R) is included and flagged.
#  NOTE: predicted PAR blends the diffuse term exactly like the pipeline via
#  eff_light_fn(); ambient is geometry-only (fixed per sensor), so the sweep
#  only moves `shadow`.
# ══════════════════════════════════════════════════════════════════════════════
buf_radii <- sort(unique(c(0.005, 0.01, 0.02, 0.05, 0.10, BUFFER_R)))  # incl. pipeline buffer

# precompute this sensor's surface_par, ambient_fraction, water_atten per hour so
# each swept shadow value turns into a correctly-blended predicted PAR.
sensor_env <- validation %>%
  filter(sitekey == pick_key, hoboID_env == pick_id) %>%
  select(datetime_pos, surface_par, ambient_fraction, water_atten, light_mol)

sweep_rows <- list()
for (k in seq_len(nrow(sun_steps))) {
  trow  <- sun_steps[k, ]
  k_alt <- trow$sun_altitude * 180 / pi
  k_az  <- trow$sun_azimuth  * 180 / pi
  
  k_shadow <- ray_shade(dem$matrix, sunaltitude = k_alt, sunangle = k_az,
                        zscale = dem$zscale, maxsearch = MAXSEARCH_PX)
  k_rast <- rast(t(k_shadow))
  ext(k_rast) <- ext(dem$rast); crs(k_rast) <- crs(dem$rast)
  
  env_k <- sensor_env[sensor_env$datetime_pos == trow$datetime_pos, ]
  
  for (br in buf_radii) {
    if (br <= 0.006) {
      # "single pixel": value at the sensor centre cell
      sv <- terra::extract(k_rast, st_coordinates(this_hobo))[1, 1]
    } else {
      bufp <- st_buffer(this_hobo, dist = br)
      sv   <- terra::extract(k_rast, vect(bufp), fun = mean, na.rm = TRUE)[, 2]
    }
    # blended effective light, matching the pipeline (ambient fixed per hour)
    pred <- if (nrow(env_k)) {
      eff <- eff_light_fn(sv, env_k$ambient_fraction[1])
      env_k$surface_par[1] * eff * env_k$water_atten[1]
    } else NA
    sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
      time      = format(trow$datetime_pos, "%m-%d %H:%M"),
      datetime  = trow$datetime_pos,
      buf_cm    = br * 100,
      is_pipe   = isTRUE(all.equal(br, BUFFER_R)),
      shadow    = sv,
      pred_par  = pred,
      measured  = if (nrow(env_k)) env_k$light_mol[1] else NA
    )
  }
}
sweep <- do.call(rbind, sweep_rows)

# ── table: shadow value by buffer radius, one row per hour ────────────────────
shadow_wide <- reshape(
  sweep[, c("time", "buf_cm", "shadow")],
  idvar = "time", timevar = "buf_cm", direction = "wide"
)
names(shadow_wide) <- gsub("shadow\\.", "buf", names(shadow_wide))
cat("\n=====================================================\n")
cat("SHADOW VALUE vs BUFFER RADIUS —", pick_id, "\n")
cat("(columns = buffer radius in cm; 0.5 cm ~ single pixel)\n")
cat("=====================================================\n")
print(shadow_wide, row.names = FALSE, digits = 3)

# ── plot: predicted PAR per buffer size vs measured, across the day ───────────
dev.new(width = 9, height = 6)
par(mar = c(7, 4, 4, 2))

times_u <- unique(sweep$datetime)
xt      <- seq_along(times_u)
cols    <- hcl.colors(length(buf_radii), "viridis", rev = TRUE)

plot(NULL, xlim = range(xt), ylim = range(c(sweep$pred_par, sweep$measured), na.rm = TRUE),
     xaxt = "n", xlab = "", ylab = "PAR (umol)",
     main = paste0(pick_id, " — predicted PAR by buffer radius vs measured"))
axis(1, at = xt, labels = format(times_u, "%m-%d %H:%M"), las = 2, cex.axis = 0.7)

# measured (thick black)
meas_line <- sweep[sweep$buf_cm == max(sweep$buf_cm), ]
lines(xt, meas_line$measured[match(times_u, meas_line$datetime)],
      col = "black", lwd = 3)
points(xt, meas_line$measured[match(times_u, meas_line$datetime)],
       col = "black", pch = 16)

# one coloured line per buffer radius
for (j in seq_along(buf_radii)) {
  br  <- buf_radii[j]
  sub <- sweep[sweep$buf_cm == br * 100, ]
  yy  <- sub$pred_par[match(times_u, sub$datetime)]
  lines(xt, yy, col = cols[j], lwd = 2, lty = 2)
  points(xt, yy, col = cols[j], pch = 1)
}

legend("topright", bty = "n", cex = 0.85,
       legend = c("Measured", paste0(buf_radii * 100, " cm buffer")),
       col    = c("black", cols),
       lwd    = c(3, rep(2, length(buf_radii))),
       lty    = c(1, rep(2, length(buf_radii))),
       pch    = c(16, rep(1, length(buf_radii))))

cat("\nBuffer sweep done. If the coloured lines sit on top of each other, buffer\n")
cat("size barely matters here. If the smallest (single-pixel) line is much\n")
cat("spikier/lower, the 5 cm buffer IS smoothing real microscale shadow.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  6. POSITION-OFFSET TEST (systematic vs random sensor mislocation)
#  For each sensor, shifts its mapped location across a grid of x/y offsets,
#  recomputes the prediction at each shifted spot, and finds the offset that
#  best matches the measured curve. Then asks: do all sensors agree on one
#  offset (systematic coordinate bug) or scatter randomly (per-sensor noise)?
#
#  How to read the output:
#   - best offsets cluster near ONE non-zero vector  -> systematic shift bug
#   - best offsets scatter in random directions      -> per-sensor position noise
#   - best offset is ~ (0,0) for most sensors         -> positions fine; look at Kd / diffuse light
# ══════════════════════════════════════════════════════════════════════════════

# --- settings -----------------------------------------------------------------
OFFSET_STEP   <- 0.05                 # grid spacing in metres (5 cm)
OFFSET_N      <- 3                    # steps each way -> grid is (2N+1)^2 = 49 offsets
offset_buf    <- 0.02                 # small buffer when sampling each shifted point (m)
offsets       <- seq(-OFFSET_N, OFFSET_N) * OFFSET_STEP   # e.g. -0.15..0.15

# sensors to test: ONLY this site's loggers (hobo_utm/dem are this site's objects,
# and sensor IDs aren't unique across sites). Cap for speed if needed.
test_ids <- sort(unique(validation$hoboID_env[validation$sitekey == pick_key]))
# test_ids <- head(test_ids, 12)      # uncomment to limit while experimenting

# precompute one shadow raster per timestep (reused across all sensors/offsets)
cat("\nPrecomputing shadow rasters for offset test (", nrow(sun_steps), " hours)...\n", sep = "")
shadow_rasters <- vector("list", nrow(sun_steps))
for (k in seq_len(nrow(sun_steps))) {
  trow  <- sun_steps[k, ]
  sh <- ray_shade(dem$matrix,
                  sunaltitude = trow$sun_altitude * 180 / pi,
                  sunangle    = trow$sun_azimuth  * 180 / pi,
                  zscale = dem$zscale, maxsearch = MAXSEARCH_PX)
  r <- rast(t(sh)); ext(r) <- ext(dem$rast); crs(r) <- crs(dem$rast)
  shadow_rasters[[k]] <- r
}

# helper: RMSE between predicted and measured for one sensor at one (dx,dy) offset
rmse_at_offset <- function(id, dx, dy) {
  base_xy <- st_coordinates(hobo_utm[hobo_utm$hoboID == id, ])
  env_s   <- validation %>%
    filter(sitekey == pick_key, hoboID_env == id) %>%
    select(datetime_pos, surface_par, ambient_fraction, water_atten, light_mol)
  if (!nrow(env_s)) return(NA_real_)
  
  pred <- meas <- numeric(0)
  for (k in seq_len(nrow(sun_steps))) {
    ek <- env_s[env_s$datetime_pos == sun_steps$datetime_pos[k], ]
    if (!nrow(ek)) next
    pt  <- cbind(base_xy[1] + dx, base_xy[2] + dy)
    sv  <- terra::extract(shadow_rasters[[k]], pt)[1, 1]   # single shifted pixel
    eff <- eff_light_fn(sv, ek$ambient_fraction[1])        # blend like pipeline
    pred <- c(pred, ek$surface_par[1] * eff * ek$water_atten[1])
    meas <- c(meas, ek$light_mol[1])
  }
  if (length(pred) < 2) return(NA_real_)
  sqrt(mean((pred - meas)^2, na.rm = TRUE))
}

# sweep every sensor over the offset grid
cat("Scanning ", length(test_ids), " sensors over a ",
    length(offsets), "x", length(offsets), " offset grid...\n", sep = "")
best <- lapply(test_ids, function(id) {
  grid <- expand.grid(dx = offsets, dy = offsets)
  grid$rmse <- mapply(function(a, b) rmse_at_offset(id, a, b), grid$dx, grid$dy)
  base_rmse <- grid$rmse[grid$dx == 0 & grid$dy == 0]
  bi <- which.min(grid$rmse)
  data.frame(hoboID = id,
             best_dx = grid$dx[bi], best_dy = grid$dy[bi],
             best_rmse = grid$rmse[bi],
             base_rmse = base_rmse,
             improve   = base_rmse - grid$rmse[bi])
})
best <- do.call(rbind, best)

cat("\n=====================================================\n")
cat("BEST POSITION OFFSET PER SENSOR\n")
cat("(best_dx/dy in metres; improve = RMSE drop vs no shift)\n")
cat("=====================================================\n")
print(best, row.names = FALSE, digits = 3)

# summary: is there a common direction?
cat("\nMean best offset across sensors:  dx =", round(mean(best$best_dx), 3),
    " dy =", round(mean(best$best_dy), 3), "m\n")
cat("SD of best offset:                dx =", round(sd(best$best_dx), 3),
    " dy =", round(sd(best$best_dy), 3), "m\n")
cat("Sensors whose best offset is (0,0):",
    sum(best$best_dx == 0 & best$best_dy == 0), "/", nrow(best), "\n")
cat("--> mean far from 0 AND small SD  = systematic shift\n")
cat("--> mean near 0 but large SD      = random per-sensor noise\n")
cat("--> most at (0,0)                 = positions fine; problem is Kd/diffuse\n")

# --- plot: arrow per sensor showing its best shift direction ------------------
dev.new(width = 7, height = 7)
lim <- max(abs(offsets)) * 1.1
plot(NULL, xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1,
     xlab = "best dx (m, +east)", ylab = "best dy (m, +north)",
     main = "Best position offset per sensor\n(clustered=systematic, scattered=random)")
abline(h = 0, v = 0, col = "grey80")
# colour by how much the shift helped
imp <- best$improve
imp_col <- hcl.colors(100, "viridis")[
  pmax(1, pmin(100, round(99 * (imp - min(imp, na.rm=T)) /
                            (diff(range(imp, na.rm=T)) + 1e-9)) + 1))]
arrows(0, 0, best$best_dx, best$best_dy, length = 0.07, col = imp_col, lwd = 2)
points(best$best_dx, best$best_dy, pch = 16, col = imp_col)
text(best$best_dx, best$best_dy, best$hoboID, pos = 3, cex = 0.6)
# mean vector in red
arrows(0, 0, mean(best$best_dx), mean(best$best_dy),
       length = 0.12, col = "red", lwd = 3)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("per-sensor best shift", "mean shift (red)"),
       col = c("grey40", "red"), lwd = c(2, 3))

cat("\nOffset test done. Look at the arrow plot: a tight fan pointing one way\n")
cat("means systematic; a starburst in all directions means random per-sensor.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  7. IS THE OFFSET A POSITION BUG OR A SHADOW-CASTING BUG?
#  Compares the mean best-offset direction (from section 6) against the
#  direction shadows are cast during the day. If the offset lines up with the
#  shadow axis, the shadows are probably mis-cast (azimuth / transpose / DEM
#  orientation) and shifting sensors just papers over it. If the offset is
#  unrelated to the sun direction, it's more likely a real coordinate bug.
#
#  Needs `best` and `sun_steps` from section 6.
# ══════════════════════════════════════════════════════════════════════════════

# --- mean offset direction (the thing section 6 found) ------------------------
mdx <- mean(best$best_dx); mdy <- mean(best$best_dy)
offset_mag <- sqrt(mdx^2 + mdy^2)

# compass bearing of the offset vector (0 = north, 90 = east), in degrees
# atan2(east, north) gives bearing clockwise from north
offset_bearing <- (atan2(mdx, mdy) * 180 / pi) %% 360

# --- shadow-cast direction --------------------------------------------------
# suncalc azimuth convention: measured from SOUTH, clockwise (so +90 = west).
# rayshader sunangle is fed az_deg = azimuth*180/pi directly in your pipeline.
# The shadow is cast in the direction OPPOSITE the sun (away from the sun).
# Convert each hour's sun azimuth to a compass bearing of the SHADOW.
sun_az_deg <- sun_steps$sun_azimuth * 180 / pi          # suncalc: from south, +CW
# bearing of the SUN (compass, from north): south=180, so sun_bearing = 180 + az
sun_bearing    <- (180 + sun_az_deg) %% 360
# shadow points away from the sun -> +180
shadow_bearing <- (sun_bearing + 180) %% 360

# weight each hour's shadow direction by how strongly it should matter:
# low sun = long shadows = bigger positional leverage. weight by cos(altitude)
# (equivalently sin of zenith) so near-horizon hours dominate, like real shadows.
alt_deg <- sun_steps$sun_altitude * 180 / pi
w       <- cos(alt_deg * pi / 180)                      # 0 overhead, ->1 near horizon
# circular weighted mean of shadow bearing
sb_rad  <- shadow_bearing * pi / 180
mean_shadow_bearing <- (atan2(sum(w * sin(sb_rad)), sum(w * cos(sb_rad))) *
                          180 / pi) %% 360

# angular difference between offset direction and shadow axis (0..180; treat as axis)
ang_diff <- abs((offset_bearing - mean_shadow_bearing + 180) %% 360 - 180)
ang_diff_axis <- min(ang_diff, 180 - ang_diff)   # shadows are an axis, not a ray

cat("\n=====================================================\n")
cat("OFFSET DIRECTION vs SHADOW DIRECTION\n")
cat("=====================================================\n")
cat("Mean best offset:           dx =", round(mdx, 3), " dy =", round(mdy, 3),
    " m  (magnitude", round(offset_mag, 3), "m)\n")
cat("Offset bearing:            ", round(offset_bearing, 1),
    "deg (compass, 0=N 90=E)\n")
cat("Mean shadow-cast bearing:  ", round(mean_shadow_bearing, 1),
    "deg (altitude-weighted)\n")
cat("Angular diff (as axis):    ", round(ang_diff_axis, 1), "deg\n\n")

if (offset_mag < OFFSET_STEP) {
  cat("VERDICT: offset magnitude <", OFFSET_STEP,
      "m (one grid step) -- essentially no net shift.\n")
  cat("         Positions look fine; chase Kd / diffuse light instead.\n")
} else if (ang_diff_axis < 30) {
  cat("VERDICT: offset is ALIGNED with the shadow axis (<30 deg).\n")
  cat("         Likely a SHADOW-CASTING problem (sun azimuth sign, matrix\n")
  cat("         transpose, or DEM orientation) -- NOT sensor misplacement.\n")
  cat("         Shifting sensors would just compensate for mis-cast shadows.\n")
} else if (ang_diff_axis > 60) {
  cat("VERDICT: offset is UNRELATED to the shadow axis (>60 deg).\n")
  cat("         More consistent with a real COORDINATE/position issue\n")
  cat("         (datum, projection params, or extent/origin mismatch).\n")
} else {
  cat("VERDICT: ambiguous (30-60 deg). Could be a mix. Inspect per-sensor\n")
  cat("         arrows in section 6 and the spread below.\n")
}

# how consistent is the offset direction across sensors? (axis-aware)
brg <- (atan2(best$best_dx, best$best_dy) * 180 / pi) %% 360
# fold to a 0..180 axis so opposite arrows count as aligned
brg_axis <- brg %% 180
# circular SD on the doubled angle (Mardia) as a simple spread measure
ca <- cos(2 * brg_axis * pi / 180); sa <- sin(2 * brg_axis * pi / 180)
R  <- sqrt(mean(ca)^2 + mean(sa)^2)         # 1 = perfectly aligned, 0 = uniform
cat("\nDirection consistency across sensors (1=aligned, 0=random):",
    round(R, 2), "\n")
cat("  (high R + alignment with shadow axis = strong shadow-casting signal;\n")
cat("   low R = scattered = random per-sensor position noise)\n")

# --- plot: offset arrows with the shadow axis overlaid ------------------------
dev.new(width = 7, height = 7)
lim <- max(abs(offsets)) * 1.1
plot(NULL, xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1,
     xlab = "dx (m, +east)", ylab = "dy (m, +north)",
     main = "Offset arrows vs shadow axis")
abline(h = 0, v = 0, col = "grey85")

# shadow axis (grey dashed line through origin, both directions)
sb <- mean_shadow_bearing * pi / 180
sx <- sin(sb) * lim; sy <- cos(sb) * lim
segments(-sx, -sy, sx, sy, col = "grey50", lwd = 2, lty = 2)
text(sx * 0.9, sy * 0.9, "shadow axis", col = "grey40", cex = 0.8)

# per-sensor offset arrows
arrows(0, 0, best$best_dx, best$best_dy, length = 0.06, col = "steelblue", lwd = 1.5)
# mean offset in red
arrows(0, 0, mdx, mdy, length = 0.12, col = "red", lwd = 3)
text(mdx, mdy, "mean offset", col = "red", pos = 3, cex = 0.8)

legend("topleft", bty = "n", cex = 0.8,
       legend = c("per-sensor offset", "mean offset", "shadow axis"),
       col = c("steelblue", "red", "grey50"),
       lwd = c(1.5, 3, 2), lty = c(1, 1, 2))

cat("\nDirection check done. If the red mean-offset arrow lies along the grey\n")
cat("shadow axis, suspect shadow-casting; if perpendicular/unrelated, suspect\n")
cat("a coordinate bug. If the arrow is tiny, positions are fine.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  8. CO-FIT Kd + AMBIENT (DIFFUSE) LIGHT
#  Position is ruled out (sections 6-7). The two remaining levers are:
#    (a) Kd too low      -> midday overshoot when shadow~1
#    (b) no diffuse floor -> deep-shadow undershoot (water scatters light in
#        from all directions, so even a fully shaded pixel sees some light)
#
#  New model:
#    eff_shadow = ambient + (1 - ambient) * shadow_fraction
#    pred_par   = surface_par * eff_shadow * exp(-Kd * depth)
#  ambient=0 recovers your current model; ambient=0.2 means shaded spots still
#  get 20% of open-sky light. Grid-searches (Kd, ambient) jointly on RMSE,
#  pooled over every sensor/hour in `validation`.
# ══════════════════════════════════════════════════════════════════════════════

# work from the already-computed shadow fractions in validation (no re-raytrace)
fitdat <- validation %>%
  filter(!is.na(shadow_fraction), !is.na(light_mol),
         !is.na(surface_par), !is.na(water_depth_m))

kd_grid  <- seq(0.05, 0.60, by = 0.01)
amb_grid <- seq(0.00, 0.40, by = 0.02)

# RMSE over the joint grid
rmse_grid <- matrix(NA_real_, length(kd_grid), length(amb_grid),
                    dimnames = list(paste0("kd", kd_grid),
                                    paste0("amb", amb_grid)))
for (a in seq_along(kd_grid)) {
  for (b in seq_along(amb_grid)) {
    eff <- amb_grid[b] + (1 - amb_grid[b]) * fitdat$shadow_fraction
    pred <- fitdat$surface_par * eff * exp(-kd_grid[a] * fitdat$water_depth_m)
    rmse_grid[a, b] <- sqrt(mean((pred - fitdat$light_mol)^2, na.rm = TRUE))
  }
}

bi   <- which(rmse_grid == min(rmse_grid), arr.ind = TRUE)[1, ]
best_kd  <- kd_grid[bi[1]]
best_amb <- amb_grid[bi[2]]

# baseline: your current model (ambient = 0) at its own best Kd, for comparison
base_col   <- rmse_grid[, which(amb_grid == 0)]
base_kd    <- kd_grid[which.min(base_col)]
base_rmse  <- min(base_col)

cat("\n=====================================================\n")
cat("JOINT Kd + AMBIENT FIT (pooled over all sensors/hours)\n")
cat("=====================================================\n")
cat("Current-style model (ambient=0):  best Kd =", base_kd,
    " RMSE =", round(base_rmse, 1), "\n")
cat("With ambient term:                best Kd =", best_kd,
    " ambient =", best_amb,
    " RMSE =", round(min(rmse_grid), 1), "\n")
cat("RMSE improvement from ambient:    ",
    round(base_rmse - min(rmse_grid), 1),
    " (", round(100 * (base_rmse - min(rmse_grid)) / base_rmse, 1), "%)\n", sep = "")

# --- heatmap of the RMSE surface ----------------------------------------------
dev.new(width = 8, height = 6)
par(mar = c(5, 5, 4, 6))
image(kd_grid, amb_grid, rmse_grid,
      col = hcl.colors(100, "viridis"),
      xlab = "Kd (m^-1)", ylab = "ambient (diffuse floor)",
      main = "RMSE surface: Kd vs ambient light")
contour(kd_grid, amb_grid, rmse_grid, add = TRUE, col = "white", nlevels = 12)
points(best_kd, best_amb, pch = 4, col = "red", cex = 2, lwd = 3)
text(best_kd, best_amb, paste0("  Kd=", best_kd, " amb=", best_amb),
     col = "red", pos = 4, cex = 0.8)

# --- before/after scatter: predicted vs measured ------------------------------
dev.new(width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

pred_old <- fitdat$surface_par * fitdat$shadow_fraction *
  exp(-base_kd * fitdat$water_depth_m)
pred_new <- fitdat$surface_par *
  (best_amb + (1 - best_amb) * fitdat$shadow_fraction) *
  exp(-best_kd * fitdat$water_depth_m)
lim <- range(c(fitdat$light_mol, pred_old, pred_new), na.rm = TRUE)

plot(fitdat$light_mol, pred_old, pch = 16, col = "#4444aa55",
     xlim = lim, ylim = lim, xlab = "Measured PAR", ylab = "Predicted PAR",
     main = paste0("Current model (amb=0, Kd=", base_kd, ")\nRMSE=",
                   round(base_rmse, 1)))
abline(0, 1, lty = 2, col = "red", lwd = 2)

plot(fitdat$light_mol, pred_new, pch = 16, col = "#aa444455",
     xlim = lim, ylim = lim, xlab = "Measured PAR", ylab = "Predicted PAR",
     main = paste0("With ambient (amb=", best_amb, ", Kd=", best_kd, ")\nRMSE=",
                   round(min(rmse_grid), 1)))
abline(0, 1, lty = 2, col = "red", lwd = 2)
par(mfrow = c(1, 1))

cat("\nIf the right-hand scatter hugs the 1:1 line tighter (esp. low-light points\n")
cat("lifting off the x-axis and high-light points pulled down), the ambient+Kd\n")
cat("combination is doing real work. Plug best_kd / best_amb into the pipeline:\n")
cat("  pred_par = surface_par * (AMBIENT + (1-AMBIENT)*shadow_fraction) *\n")
cat("             exp(-KD_CONSTANT * water_depth_m)\n")


# ══════════════════════════════════════════════════════════════════════════════
#  9. WHY IS THE MODEL WORST AT HIGH SUN?
#  Hypothesis: when the sun is near-vertical, shadows collapse, every pixel
#  returns ~1, the shadow term stops separating sensors, AND ray_shade becomes
#  hypersensitive to zscale. This section tests both.
# ══════════════════════════════════════════════════════════════════════════════

# ── 9a. error vs sun altitude — is midday really the worst? ───────────────────
err <- validation %>%
  filter(!is.na(pred_par), !is.na(light_mol)) %>%
  mutate(alt_deg  = sun_altitude * 180 / pi,
         resid    = pred_par - light_mol,
         abs_resid = abs(resid))

# bin by sun altitude and summarise error + how flat the shadow field is
alt_bins <- err %>%
  mutate(alt_bin = cut(alt_deg, breaks = seq(0, 90, by = 10))) %>%
  group_by(alt_bin) %>%
  summarise(
    n            = n(),
    mean_shadow  = round(mean(shadow_fraction, na.rm = TRUE), 3),
    sd_shadow    = round(sd(shadow_fraction,   na.rm = TRUE), 3),  # spread = info content
    bias         = round(mean(resid, na.rm = TRUE), 1),
    rmse         = round(sqrt(mean(resid^2, na.rm = TRUE)), 1),
    .groups = "drop"
  )
cat("\n=====================================================\n")
cat("ERROR vs SUN ALTITUDE\n")
cat("(sd_shadow = how much shadow varies between sensors;\n")
cat(" low sd at high sun = shadow term has stopped separating loggers)\n")
cat("=====================================================\n")
print(as.data.frame(alt_bins), row.names = FALSE)

dev.new(width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
plot(err$alt_deg, err$resid, pch = 16, col = "#33333333",
     xlab = "Sun altitude (deg)", ylab = "Predicted - Measured PAR",
     main = "Residual vs sun altitude")
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(lowess(err$alt_deg, err$resid, f = 0.5), col = "blue", lwd = 2)

# shadow spread vs altitude: does the shadow field go flat at high sun?
plot(err$alt_deg, err$shadow_fraction, pch = 16, col = "#33aa3333",
     xlab = "Sun altitude (deg)", ylab = "shadow_fraction",
     main = "Shadow value vs sun altitude\n(collapse toward 1 = losing info)")
lines(lowess(err$alt_deg, err$shadow_fraction, f = 0.5), col = "darkgreen", lwd = 2)
par(mfrow = c(1, 1))


# ── 9b. zscale sensitivity at high vs low sun ─────────────────────────────────
# Re-raytrace ONE high-sun hour and ONE low-sun hour across a range of zscales,
# and watch how much each sensor's shadow value moves. If high sun is far more
# sensitive, zscale error is contaminating your midday predictions specifically.
zscale_mult <- c(0.5, 0.75, 1, 1.5, 2, 3)        # multiples of current dem$zscale

hi <- sun_steps[which.max(sun_steps$sun_altitude), ]
lo <- sun_steps[which.min(sun_steps$sun_altitude), ]

probe_zscale <- function(trow) {
  alt <- trow$sun_altitude * 180 / pi
  az  <- trow$sun_azimuth  * 180 / pi
  sapply(zscale_mult, function(m) {
    sh <- ray_shade(dem$matrix, sunaltitude = alt, sunangle = az,
                    zscale = dem$zscale * m, maxsearch = MAXSEARCH_PX)
    r <- rast(t(sh)); ext(r) <- ext(dem$rast); crs(r) <- crs(dem$rast)
    # mean shadow across all sensor buffers
    bv <- terra::extract(r, vect(st_buffer(hobo_utm, BUFFER_R)), fun = mean, na.rm = TRUE)[, 2]
    mean(bv, na.rm = TRUE)
  })
}

hi_curve <- probe_zscale(hi)
lo_curve <- probe_zscale(lo)

cat("\n=====================================================\n")
cat("ZSCALE SENSITIVITY (mean sensor shadow vs zscale multiplier)\n")
cat("current zscale =", round(dem$zscale, 5), "m\n")
cat("=====================================================\n")
zs_tbl <- data.frame(zscale_mult = zscale_mult,
                     high_sun = round(hi_curve, 3),
                     low_sun  = round(lo_curve, 3))
print(zs_tbl, row.names = FALSE)
cat("\nhigh_sun spread:", round(diff(range(hi_curve)), 3),
    " | low_sun spread:", round(diff(range(lo_curve)), 3), "\n")
cat("If high_sun moves much more across zscales than low_sun, midday shadows\n")
cat("are zscale-sensitive and a wrong zscale is poisoning your noon predictions.\n")

dev.new(width = 7, height = 6)
par(mar = c(5, 5, 4, 2))
plot(zscale_mult, hi_curve, type = "b", pch = 16, col = "orange", lwd = 2,
     ylim = range(c(hi_curve, lo_curve)),
     xlab = "zscale multiplier (1 = current)", ylab = "mean sensor shadow fraction",
     main = paste0("zscale sensitivity: high sun (",
                   round(hi$sun_altitude*180/pi), "deg) vs low sun (",
                   round(lo$sun_altitude*180/pi), "deg)"))
lines(zscale_mult, lo_curve, type = "b", pch = 17, col = "navy", lwd = 2)
abline(v = 1, lty = 3, col = "grey50")
legend("topright", c("high sun (noon)", "low sun (morning/eve)"),
       col = c("orange", "navy"), pch = c(16, 17), lwd = 2, bty = "n")

cat("\nHigh-sun diagnosis done.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  10. DOWNSAMPLE SENSITIVITY — does fine resolution recover midday shadow?
#  At fact=10, ~1mm native pixels become ~1cm. Fine rugosity (coral branches,
#  small overhangs) that casts real midday shadows gets averaged away, so the
#  model CAN'T produce noon variation even where it physically exists.
#  This rebuilds the DEM at several downsample factors and, for one high-sun
#  hour, measures how much sensor-to-sensor shadow VARIATION (sd) returns as
#  resolution improves. Rising sd at finer resolution = real midday structure
#  was being smoothed out.
#
#  WARNING: small fact = big matrix = slow + memory-hungry. Start conservative.
#  Uses the picked sensor's site/year (pick_key). Single site/year only.
#
#  ⚠ API NOTE: the current master `prepare_dem(site, year, target_px)` has no
#  `fact =` argument — downsampling is driven by TARGET_PIXELS and the factor is
#  derived per DEM. We therefore convert each desired `fact` into the equivalent
#  target_px (= native_px / fact) and pass that. A true fixed-factor sweep would
#  need a fact-based prepare_dem variant.
# ══════════════════════════════════════════════════════════════════════════════

# resolve the picked site/year from pick_key via dem_store
active_site <- dem_store[[pick_key]]$site
active_year <- dem_store[[pick_key]]$year

# guard: this needs a single active site/year (prepare_dem takes one each)
if (length(active_site) != 1 || length(active_year) != 1) {
  cat("\n[Section 10 skipped] pick_key did not resolve to one site/year.\n")
} else {
  
  # factors to test: smaller = finer = slower.
  # edit downward (e.g. add 3, 2) once you know your machine can handle it.
  ds_factors <- c(20, 10, 5)
  
  # one high-sun hour (worst case for the model, per section 9)
  hi_row <- sun_steps[which.max(sun_steps$sun_altitude), ]
  hi_alt <- hi_row$sun_altitude * 180 / pi
  hi_az  <- hi_row$sun_azimuth  * 180 / pi
  
  cat("\n=====================================================\n")
  cat("DOWNSAMPLE SENSITIVITY at high sun (",
      round(hi_alt), "deg ) —", active_site, active_year, "\n")
  cat("=====================================================\n")
  
  # native (pre-downsample) pixel count, to turn `fact` into a target_px.
  native_px_full <- max(nrow(rast(file.path(DEM_DIR, paste0(active_site, active_year, ".tif")))),
                        ncol(rast(file.path(DEM_DIR, paste0(active_site, active_year, ".tif")))))
  
  ds_rows <- list()
  for (f in ds_factors) {
    cat("  building DEM at fact ~", f, "...\n")
    dchk <- tryCatch(
      prepare_dem(active_site, active_year,
                  reef_radius_px = max(50, round(REEF_RADIUS_PX_TARGET * (10 / f)))),
      error = function(e) { cat("    failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(dchk)) next
    
    sh <- ray_shade(dchk$matrix, sunaltitude = hi_alt, sunangle = hi_az,
                    zscale = dchk$zscale, maxsearch = MAXSEARCH_PX)
    r <- rast(t(sh)); ext(r) <- ext(dchk$rast); crs(r) <- crs(dchk$rast)
    
    # sample all sensor buffers; record spread (sd) and mean
    bv <- terra::extract(r, vect(st_buffer(hobo_utm, BUFFER_R)),
                         fun = mean, na.rm = TRUE)[, 2]
    ds_rows[[length(ds_rows) + 1]] <- data.frame(
      fact         = f,
      pixel_m      = round(dchk$zscale, 4),     # zscale == pixel size here
      matrix_px    = nrow(dchk$matrix),
      mean_shadow  = round(mean(bv, na.rm = TRUE), 3),
      sd_shadow    = round(sd(bv,   na.rm = TRUE), 3),
      frac_shaded  = round(mean(bv < 0.5, na.rm = TRUE), 3)  # share of sensors in real shade
    )
  }
  ds_tbl <- do.call(rbind, ds_rows)
  ds_tbl <- ds_tbl[order(ds_tbl$fact, decreasing = TRUE), ]
  print(ds_tbl, row.names = FALSE)
  
  cat("\nsd_shadow = sensor-to-sensor variation at this hour.\n")
  cat("If sd_shadow RISES as fact drops (finer resolution), then midday shadow\n")
  cat("structure IS being smoothed away at fact=10 — finer DEM would help.\n")
  cat("If sd_shadow stays flat, resolution isn't the midday limiter; the reef\n")
  cat("genuinely casts little shadow at noon and the model is at its honest floor.\n")
  
  # plot
  dev.new(width = 8, height = 6)
  par(mar = c(5, 5, 4, 5))
  o <- order(ds_tbl$pixel_m)
  plot(ds_tbl$pixel_m[o], ds_tbl$sd_shadow[o], type = "b", pch = 16,
       col = "firebrick", lwd = 2,
       xlab = "pixel size (m)  [<- finer    coarser ->]",
       ylab = "sd of sensor shadow values",
       main = paste0("Does midday shadow variation return at finer resolution?\n",
                     active_site, active_year, " @ ", round(hi_alt), " deg sun"))
  # secondary: mean shadow on right axis
  par(new = TRUE)
  plot(ds_tbl$pixel_m[o], ds_tbl$mean_shadow[o], type = "b", pch = 1, lty = 2,
       col = "grey50", axes = FALSE, xlab = "", ylab = "")
  axis(4); mtext("mean shadow (dashed)", side = 4, line = 3, cex = 0.9)
  abline(v = dem$zscale, lty = 3, col = "blue")
  text(dem$zscale, min(ds_tbl$mean_shadow),
       paste0("current\n(", round(dem$zscale, 3), " m)"),
       col = "blue", cex = 0.7, pos = 4)
  
  cat("\nDownsample test done. Restore dem/hobo to your chosen fact afterward if\n")
  cat("you changed anything (this section uses local DEMs, doesn't touch `dem`).\n")
}