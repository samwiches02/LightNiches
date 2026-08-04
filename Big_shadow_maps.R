# ══════════════════════════════════════════════════════════════════════════════
#  BIG SHADOW MAPS — one PNG per timestep: shadow map on top, pred vs measured below
# ══════════════════════════════════════════════════════════════════════════════
sun_steps <- output %>%
  distinct(datetime_pos, sun_altitude, sun_azimuth) %>%
  arrange(datetime_pos)

out_dir <- file.path(getwd(), "shadow_frames")
dir.create(out_dir, showWarnings = FALSE)

site_tag <- paste0(unique(output$site[!is.na(output$site)]),
                   unique(output$year[!is.na(output$year)]))

# shared y-limit for the bar charts so frames are comparable
bar_ymax <- max(c(validation$pred_par, validation$light_mol), na.rm = TRUE)

for (i in seq_len(nrow(sun_steps))) {
  ts      <- sun_steps[i, ]
  alt_deg <- ts$sun_altitude * 180 / pi
  az_deg  <- ts$sun_azimuth  * 180 / pi
  
  # raytrace this timestep (same call as the pipeline)
  shadow <- ray_shade(dem$matrix,
                      sunaltitude = alt_deg,
                      sunangle    = az_deg,
                      zscale      = dem$zscale,
                      maxsearch   = MAXSEARCH_FACT)
  shadow[!dem$reef_mask] <- NA          # blank off-reef border
  
  shadow_rast <- rast(shadow)
  ext(shadow_rast) <- ext(dem$rast)
  crs(shadow_rast) <- crs(dem$rast)     # UTM 55S
  shadow_wgs <- project(shadow_rast, "EPSG:4326")
  
  # one tall file per timestep
  fname <- file.path(
    out_dir,
    paste0(site_tag, "_", format(ts$datetime_pos, "%Y%m%d_%H%M"), ".png")
  )
  png(fname, width = 1400, height = 1900, res = 150)
  
  # two stacked panels: map (tall) on top, bar chart (short) below
  layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1.6))
  
  # ── panel 1: shadow map ──
  par(mar = c(4, 4, 5, 6))
  plot(shadow_wgs,
       main = paste0(site_tag, "  |  ",
                     format(ts$datetime_pos, "%Y-%m-%d %H:%M"),
                     "\nsun alt = ", round(alt_deg, 1), "°   az = ", round(az_deg, 1),
                     "°   |  maxsearch=", MAXSEARCH_FACT,
                     "  downsample=", DOWNSAMPLE_FACTOR),
       cex.main = 1.0)
  
  plot(vect(hobo_wgs84), add = TRUE,
       col = "red", bg = NA, cex = 1.8, pch = 21, lwd = 2)
  text(st_coordinates(hobo_wgs84)[, 1],
       st_coordinates(hobo_wgs84)[, 2],
       labels = hobo_sf$hoboID,
       col = "yellow", cex = 1.0, font = 2, pos = 3)
  
  # ── panel 2: predicted vs measured bar chart for this hour ──
  par(mar = c(5, 4, 2, 6))
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
            cex.names = 0.7,
            las       = 2,
            ylim      = c(0, bar_ymax),
            border    = NA,
            main      = "Predicted (red) vs Measured (blue) PAR",
            cex.main  = 0.9)
    legend("topright",
           legend = c("Predicted", "Measured"),
           fill   = c("red", "blue"), border = NA, bty = "n", cex = 0.9)
  } else {
    plot.new()
    text(0.5, 0.5, "no validation data this hour", cex = 1.1)
  }
  
  dev.off()
  message("saved ", basename(fname))
}

layout(1)
message("All frames in: ", out_dir)
