# ══════════════════════════════════════════════════════════════════════════════
#  04_figures.R — figures for the paper and supplement
#
#  Reads : output/predictions.rds, output/metrics_*.csv, output/ref_number.csv,
#          output/cache.rds, output/dems/<FIG_SITE>.rds, output/geom.rds
#  Writes: output/figures/*.png
#
#  Figure numbers follow the L&O: Methods skeleton:
#    Fig 3  shadow maps, in-air vs refracted sun        (needs 02 outputs only)
#    Fig 4  model ladder: RMSE and bias with 95% CIs
#    Fig 5  within-hour anomalies, measured vs predicted, by event
#    Fig 6  error vs number of reference loggers, open vs random
#    Fig 7  residuals vs hour, depth and shadow fraction
#    S1     per-logger time series per event
#    S6     daytime-integrated PAR, measured vs predicted
#
#  PROVENANCE TAGS:  [YOURS] your code from master_script.R, logic unchanged
#                    [YOURS, EDITED] your code with listed changes   [NEW] Claude
#  Comments with no "[Claude]" marker are yours, copied from master_script.R.
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_raytrace.R"))
source(here::here("R", "functions_models.R"))
suppressPackageStartupMessages({ library(dplyr) })

stage_banner("04  FIGURES")
init_geometry_checks()   # [NEW] same orientation checks and sun calibration as stage 02
fig <- function(name, w = 8, h = 6) {
  png(out_path("figures", paste0(name, ".png")), width = w, height = h, units = "in", res = 200)
}
UNIT_HOURLY <- expression("PAR (" * mu * "mol photons m"^-2 * " s"^-1 * ")")

# ── Fig 3: shadow maps, in-air vs refracted, at FIG_HOURS ────────────────────
# [YOURS, EDITED] — based on your plot_shadow_3hr_with_zoom(): same day choice
#   (day with the most timesteps), nearest timestep to each target hour,
#   grayscale palette, sensors overlaid
#   - [Claude] two rows (unrefracted vs refracted) instead of the zoom panels
#     and the full-day line; your zoom panels are not reproduced
#   - [Claude] uses the corrected angles and orientation
#  Row 1 (a): three grayscale shadow maps (08, 12, 16h), sensors overlaid
#  Grayscale = ready for print; white = lit, dark = shaded.
#  Needs only 02 outputs, so it can be drawn before the models are run.
local({
  dem   <- load_dem(out_path("dems", paste0(FIG_SITE, ".rds")))
  hobo  <- load_hobo_locations(readRDS(out_path("geom.rds")), FIG_SITE)
  cache <- readRDS(out_path("cache.rds"))
  ts    <- sort(unique(cache$datetime_pos[cache$sitekey == FIG_SITE]))
  day   <- names(sort(table(as.Date(ts, tz = LOCAL_TZ)), decreasing = TRUE))[1]
  ts    <- ts[as.Date(ts, tz = LOCAL_TZ) == as.Date(day)]
  hrs   <- as.integer(format(ts, "%H"))
  pick  <- ts[vapply(FIG_HOURS, function(h) which.min(abs(hrs - h)), integer(1))]
  cen   <- get_site_centroid(dem$rast)
  pal   <- gray.colors(100, start = 0.1, end = 0.98)

  fig("fig3_shadow_maps_refraction", w = 3.2 * length(pick), h = 7)
  par(mfrow = c(2, length(pick)), mar = c(0.5, 0.5, 3, 0.5), oma = c(0, 3, 2, 0))
  for (refract in c(FALSE, TRUE)) for (i in seq_along(pick)) {
    sun <- sun_at_midhour(pick[i], cen$lat, cen$lon)
    a   <- raytrace_angles(sun, refract)
    r   <- terra::trim(matrix_to_rast(make_shadow(dem, sun, refract), dem$rast))
    terra::plot(r, col = pal, legend = FALSE, axes = FALSE,
                main = sprintf("%s  sun %.0f°%s", format(pick[i] + 1800, "%H:%M"),
                               a$alt_deg, if (refract) " (in water)" else " (in air)"))
    terra::plot(terra::vect(hobo), add = TRUE, pch = 21, col = "red", cex = 0.9)
    if (i == 1) mtext(if (refract) "Refracted" else "Unrefracted", side = 2, line = 1)
  }
  mtext(sprintf("%s, %s: direct-light fraction (white = lit)", FIG_SITE, day),
        outer = TRUE, line = 0.3)
  dev.off()
})

if (!file.exists(out_path("predictions.rds"))) {
  end_stage("No predictions yet (03 not finished); only Fig 3 drawn.", 0L)
}

preds <- readRDS(out_path("predictions.rds"))
cache <- readRDS(out_path("cache.rds"))
boot  <- read.csv(out_path("metrics_bootstrap_ci.csv"), stringsAsFactors = FALSE)
refn  <- read.csv(out_path("ref_number.csv"), stringsAsFactors = FALSE)
main  <- preds %>% filter(model == FIG_MAIN_MODEL) %>%
  left_join(cache %>% select(sitekey, hoboID_env, datetime_pos, shadow_refr),
            by = c("sitekey", "hoboID_env", "datetime_pos"))
events <- sort(unique(preds$sitekey))
ev_col <- setNames(hcl.colors(length(events), "Dark 3"), events)

# ── Fig 4: ladder ────────────────────────────────────────────────────────────
# [NEW] — same idea as your thesis Fig 7, with confidence intervals
local({
  ladder <- c("M0", "M1", "M2s", "M3s", "M4s", "M2r", "M3r", "M4r")
  b <- boot[boot$model %in% ladder, ]
  fig("fig4_model_ladder", w = 9, h = 5)
  par(mfrow = c(1, 2), mar = c(4.5, 6, 2, 1))
  for (met in c("rmse", "bias")) {
    d <- b[b$metric == met, ]; d <- d[match(rev(ladder), d$model), ]
    y <- seq_len(nrow(d))
    plot(d$mid, y, xlim = range(c(d$lo, d$hi, 0)), yaxt = "n", pch = 19,
         xlab = if (met == "rmse") "Hourly RMSE" else "Mean bias (pred − meas)", ylab = "",
         main = toupper(met))
    segments(d$lo, y, d$hi, y, lwd = 2)
    axis(2, at = y, labels = d$model, las = 1)
    if (met == "bias") abline(v = 0, lty = 2)
    abline(h = 3.5, col = "grey80")   # separates satellite (s) and reference (r) branches
  }
  dev.off()
})

# ── Fig 5: within-hour anomalies ─────────────────────────────────────────────
# [NEW]
local({
  a <- within_anomalies(main)
  fig("fig5_within_site_anomalies", w = 7, h = 7)
  lim <- range(c(a$d_obs, a$d_pred))
  plot(a$d_obs, a$d_pred, pch = 16, cex = 0.5, col = adjustcolor(ev_col[a$sitekey], 0.5),
       xlim = lim, ylim = lim, xlab = "Measured − site-hour mean",
       ylab = "Predicted − site-hour mean",
       main = sprintf("%s: logger-to-logger differences (within-site R² = %.2f)",
                      FIG_MAIN_MODEL, within_r2(main)))
  abline(0, 1, lty = 2); abline(h = 0, v = 0, col = "grey80")
  legend("topleft", legend = events, col = ev_col, pch = 16, bty = "n", cex = 0.8)
  dev.off()
})

# ── Fig 6: reference-logger number ───────────────────────────────────────────
# [NEW]
local({
  s <- refn %>% group_by(mode, k) %>%
    summarise(med = stats::median(rmse), lo = stats::quantile(rmse, 0.05),
              hi = stats::quantile(rmse, 0.95), .groups = "drop")
  fig("fig6_reference_logger_number", w = 6, h = 5)
  plot(NA, xlim = range(s$k) + c(-0.3, 0.3), ylim = range(c(s$lo, s$hi)),
       xlab = "Reference loggers per site", ylab = "Hourly RMSE (M3, reference Kd)",
       xaxt = "n")
  axis(1, at = REF_KS)
  cols <- c(open = "#0072B2", random = "#D55E00")
  for (m in names(cols)) {
    d <- s[s$mode == m, ]; off <- if (m == "open") -0.08 else 0.08
    segments(d$k + off, d$lo, d$k + off, d$hi, col = cols[m], lwd = 2)
    lines(d$k + off, d$med, type = "b", pch = 19, col = cols[m])
  }
  legend("topright", legend = c("Open positions", "Random positions"),
         col = cols, pch = 19, lty = 1, bty = "n")
  dev.off()
})

# ── Fig 7: residuals ─────────────────────────────────────────────────────────
# [YOURS, EDITED] — same three panels as your thesis Fig 11 (hour, depth,
#   shadow fraction), pooled across events instead of one row per event
local({
  main$res  <- main$pred - main$light_mol
  main$hour <- as.integer(format(main$datetime_pos, "%H"))
  fig("fig7_residuals", w = 12, h = 4.5)
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1))
  boxplot(res ~ hour, data = main, xlab = "Hour (local)", ylab = "Residual (pred − meas)",
          outline = FALSE, col = "grey90"); abline(h = 0, lty = 2, col = "red")
  for (v in c("water_depth_m", "shadow_refr")) {
    plot(main[[v]], main$res, pch = 16, cex = 0.4, col = adjustcolor(ev_col[main$sitekey], 0.4),
         xlab = if (v == "water_depth_m") "Water depth (m)" else "Direct-light fraction",
         ylab = "Residual (pred − meas)")
    abline(h = 0, lty = 2, col = "red")
    ok <- is.finite(main[[v]]); lines(stats::lowess(main[[v]][ok], main$res[ok]), lwd = 2)
  }
  dev.off()
})

# ── S1: per-logger time series, one file per event ───────────────────────────
# [YOURS, EDITED] — your plot_par_panels_one_site(): one panel per logger,
#   5 columns, shared y axis, measured solid blue, predicted dashed red
#   - [Claude] title strip and legend row simplified
#  Useful to spot specific loggers where the prediction badly diverges from the
#  measurement — they usually have one of: wrong placement (lon/lat off),
#  wrong dpHB, or unmodelled overhang (3-D cavity the top-down DEM flattens).
for (k in events) {
  d <- main[main$sitekey == k, ]
  ids <- sort(unique(d$hoboID_env)); nc <- 5; nr <- ceiling(length(ids) / nc)
  fig(paste0("S1_timeseries_", k), w = 2.4 * nc, h = 2 * nr + 0.6)
  par(mfrow = c(nr, nc), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 2, 0))
  yl <- range(c(d$pred, d$light_mol))
  for (id in ids) {
    s <- d[d$hoboID_env == id, ]; s <- s[order(s$datetime_pos), ]
    plot(s$datetime_pos, s$light_mol, type = "l", col = "steelblue", lwd = 1.5,
         ylim = yl, main = id, cex.main = 0.8, xlab = "", ylab = "")
    lines(s$datetime_pos, s$pred, col = "firebrick", lty = 2, lwd = 1.5)
  }
  mtext(sprintf("%s  %s  (blue measured, red predicted)", k, FIG_MAIN_MODEL),
        outer = TRUE, line = 0.5)
  dev.off()
}

# ── S6: daytime-integrated PAR ───────────────────────────────────────────────
# [YOURS, EDITED] — your POOLED FIGURE 1 (DLI predicted vs measured)
#   - [Claude] colour = event (you had shape = event, colour = depth); sensor
#     labels dropped
local({
  di <- daytime_integrated(main)
  fig("S6_daytime_integrated", w = 6.5, h = 6)
  lim <- range(c(di$obs, di$pred))
  plot(di$obs, di$pred, pch = 19, col = ev_col[di$sitekey], xlim = lim, ylim = lim,
       xlab = expression("Measured daytime-integrated PAR (mol m"^-2 * ")"),
       ylab = expression("Predicted (mol m"^-2 * ")"))
  abline(0, 1, lty = 2)
  legend("topleft", legend = events, col = ev_col, pch = 19, bty = "n", cex = 0.8)
  dev.off()
})

write_manifest("04_figures")
message("04 done. Figures in ", out_path("figures"))
