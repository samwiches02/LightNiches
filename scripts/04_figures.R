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
#    S1     per-logger time series per event   -> figures/timeseries/
#    S6     daytime-integrated PAR, measured vs predicted
#    Fig 8  hourly PAR: baseline M1 vs primary M4r
#    Fig 9  daytime-integrated PAR: the same comparison
#    Fig 10 event curves for the best, median and worst event (M1 vs M4r)
#    S8     the same for every event               -> figures/event_curves/
#    S7     residuals by event vs hour, depth and direct-light fraction
#                                               -> figures/residuals_by_event/
#  Single-logger shadow checks: scripts/08_logger_check.R (run on demand)
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
# [Claude, 22 Oct 2026] `name` may include a subfolder ("timeseries/S1_..."), and
#   every figure gets a small grey footer: run label | run id | git commit
#   (* = uncommitted changes), so a PNG always says which run drew it.
FIG_STAMP <- sprintf("%s | run %s | commit %s%s", RUN_LABEL, RUN_ID, git_hash(),
                     if (git_dirty()) "*" else "")
.fig_open <- FALSE
fig <- function(name, w = 8, h = 6) {
  png(out_path("figures", paste0(name, ".png")), width = w, height = h, units = "in", res = 200)
  .fig_open <<- TRUE
}
dev.off <- function(which = grDevices::dev.cur()) {
  if (.fig_open) {
    graphics::mtext(FIG_STAMP, side = 1, line = -1.1, outer = TRUE, adj = 0.995,
                    cex = 0.55, col = "grey55")
    .fig_open <<- FALSE
  }
  grDevices::dev.off(which)
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
  dem   <- load_dem(in_path("dems", paste0(FIG_SITE, ".rds")))
  hobo  <- load_hobo_locations(readRDS(in_path("geom.rds")), FIG_SITE)
  cache <- readRDS(in_path("cache.rds"))
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

if (!file.exists(in_path("predictions.rds"))) {
  end_stage("No predictions yet (03 not finished); only Fig 3 drawn.", 0L)
}

preds <- readRDS(in_path("predictions.rds"))
cache <- readRDS(in_path("cache.rds"))
boot  <- read.csv(in_path("metrics_bootstrap_ci.csv"), stringsAsFactors = FALSE)
refn  <- read.csv(in_path("ref_number.csv"), stringsAsFactors = FALSE)
main  <- preds %>% filter(model == FIG_MAIN_MODEL) %>%
  left_join(cache %>% select(sitekey, hoboID_env, datetime_pos, shadow_refr),
            by = c("sitekey", "hoboID_env", "datetime_pos"))
events <- sort(unique(preds$sitekey))
ev_col <- setNames(hcl.colors(length(events), "Dark 3"), events)

# ── Fig 4: ladder ────────────────────────────────────────────────────────────
# [NEW] — same idea as your thesis Fig 7, with confidence intervals
#   [Claude, 22 Oct 2026] third panel: within-site R² (logger-to-logger skill);
#   primary model in red; satellite branch labelled with its 12 events
local({
  ladder <- c("M0", "M1", "M1d", "M2s", "M3s", "M4s", "M2r", "M3r", "M4r")
  labs <- c(M0 = "M0  surface PAR", M1 = "M1  + constant Kd", M1d = "M1d  + direct shading (baseline)",
            M2s = "M2s  satellite Kd", M3s = "M3s  + direct shading", M4s = "M4s  + ambient shading",
            M2r = "M2r  reference Kd", M3r = "M3r  + direct shading", M4r = "M4r  + ambient shading")
  b <- boot[boot$model %in% ladder, ]
  fig("fig4_model_ladder", w = 12, h = 5)
  par(mfrow = c(1, 3), mar = c(4.5, 1, 3, 1), oma = c(0, 11, 0, 0))
  panels <- c(rmse = "Hourly RMSE", bias = "Mean bias (pred - meas)",
              within_r2 = "Within-site R\u00B2")
  for (met in names(panels)) {
    d <- b[b$metric == met, ]; d <- d[match(rev(ladder), d$model), ]
    y <- seq_len(nrow(d))
    col <- ifelse(d$model == FIG_MAIN_MODEL, "#C8372D", ifelse(grepl("s$", d$model), "grey55", "black"))
    plot(d$mid, y, xlim = range(c(d$lo, d$hi, 0), na.rm = TRUE), ylim = c(0.5, length(y) + 0.5),
         yaxt = "n", pch = 19, col = col, cex = 1.3, xlab = panels[[met]], ylab = "",
         main = panels[[met]])
    segments(d$lo, y, d$hi, y, lwd = 2.5, col = col)
    if (met == names(panels)[1]) axis(2, at = y, labels = labs[d$model], las = 1, tick = FALSE,
                                      cex.axis = 1.05, outer = FALSE, line = 0, hadj = 1)
    if (met != "rmse") abline(v = 0, lty = 2, col = "grey40")
    abline(h = 3.5, col = "grey80")   # separates satellite (s) and reference (r) branches
    if (met == "rmse") {
      text(par("usr")[2], 5, "satellite Kd\n(12 events)", adj = c(1.05, 0.5), cex = 0.85, col = "grey45")
      text(par("usr")[2], 2, "reference Kd\n(17 events)", adj = c(1.05, 0.5), cex = 0.85)
    }
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
       xlim = lim, ylim = lim, xlab = "Measured - site-hour mean",
       ylab = "Predicted - site-hour mean",
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
  boxplot(res ~ hour, data = main, xlab = "Hour (local)", ylab = "Residual (pred - meas)",
          outline = FALSE, col = "grey90"); abline(h = 0, lty = 2, col = "red")
  for (v in c("water_depth_m", "shadow_refr")) {
    plot(main[[v]], main$res, pch = 16, cex = 0.4, col = adjustcolor(ev_col[main$sitekey], 0.4),
         xlab = if (v == "water_depth_m") "Water depth (m)" else "Direct-light fraction",
         ylab = "Residual (pred - meas)")
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
  fig(paste0("timeseries/S1_timeseries_", k), w = 2.4 * nc, h = 2 * nr + 0.6)
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

# ══════════════════════════════════════════════════════════════════════════════
#  Baseline (M1) vs primary model (M4r), all events pooled
#  [NEW] Claude, 22 Oct 2026
#    Baseline = M1d: constant Kd (KD_CONSTANT) + direct shading. [29 Oct 2026; was M1]
#    Primary  = M4r: reference-logger hourly Kd + direct and ambient shading.
#    fig8   hourly PAR, measured vs predicted, coloured by hour
#    fig9   daytime-integrated PAR, shape = event, colour = depth
#    S7     residuals by event: vs hour, depth, direct-light fraction
# ══════════════════════════════════════════════════════════════════════════════
cmp <- local({
  m1  <- preds %>% filter(model == "M1d") %>%
    select(sitekey, hoboID_env, datetime_pos, pred_m1 = pred)
  m4  <- preds %>% filter(model == FIG_MAIN_MODEL) %>%
    select(sitekey, hoboID_env, datetime_pos, light_mol, water_depth_m,
           pred_ref = pred, w_fold)
  out <- m4 %>%
    inner_join(m1, by = c("sitekey", "hoboID_env", "datetime_pos")) %>%
    left_join(cache %>% select(sitekey, hoboID_env, datetime_pos, shadow_refr, ambient),
              by = c("sitekey", "hoboID_env", "datetime_pos")) %>%
    mutate(pred_base = pred_m1,                                  # baseline = M1d
           hour       = as.integer(format(datetime_pos, "%H")))
  stopifnot(nrow(out) == nrow(m4), all(is.finite(out$pred_base)))
  out
})
CMP_MODELS <- c(pred_base = sprintf("Baseline (M1d): constant Kd = %.1f + direct shading", KD_CONSTANT),
                pred_ref  = sprintf("Primary (%s): reference Kd + shading", FIG_MAIN_MODEL))
stat_box <- function(obs, pred, extra = NULL) {
  m <- metrics(obs, pred)
  legend("topleft", bty = "n", cex = 0.9, legend = c(
    sprintf("n = %d", m$n), sprintf("RMSE = %.0f", m$rmse),
    sprintf("bias = %+.0f", m$bias), sprintf("R² = %.2f", m$r2), extra))
}
fmt_stat_box <- function(obs, pred, digits = 2) {   # DLI version: more decimals
  m <- metrics(obs, pred); f <- stats::coef(stats::lm(obs ~ pred))
  legend("topleft", bty = "n", cex = 0.9, legend = c(
    sprintf("n = %d", m$n), sprintf("RMSE = %.*f", digits, m$rmse),
    sprintf("bias = %+.*f", digits, m$bias), sprintf("R² = %.2f", m$r2),
    sprintf("fit: y = %.2f + %.2fx", f[1], f[2])))
}

# ── Fig 8: hourly PAR, measured vs predicted ─────────────────────────────────
local({
  hrs   <- sort(unique(cmp$hour))
  hcol  <- setNames(colorRampPalette(c("#3B0F0A", "#C8553D", "#F4C9B8",
                                       "#C8553D", "#3B0F0A"))(length(hrs)), hrs)
  show  <- intersect(c(8, 12, 16), hrs)
  lim   <- c(0, max(c(cmp$light_mol, cmp$pred_base, cmp$pred_ref)))
  fig("fig8_hourly_M1_vs_M4r", w = 12, h = 6)
  layout(matrix(c(1, 2, 3), 1), widths = c(1, 1, 0.18))
  par(mar = c(4.5, 4.8, 3, 1))
  for (p in names(CMP_MODELS)) {
    x <- cmp[[p]]; y <- cmp$light_mol
    plot(x, y, pch = 16, cex = 0.45, col = adjustcolor(hcol[as.character(cmp$hour)], 0.7),
         xlim = lim, ylim = lim, xlab = expression("Predicted PAR (" * mu * "mol m"^-2 * " s"^-1 * ")"),
         ylab = expression("Measured PAR (" * mu * "mol m"^-2 * " s"^-1 * ")"),
         main = CMP_MODELS[[p]], cex.main = 1)
    abline(0, 1, lty = 3, col = "grey40")
    sl <- vapply(show, function(h) { i <- cmp$hour == h; sum(x[i] * y[i]) / sum(x[i]^2) }, 0)
    for (j in seq_along(show)) {
      i <- cmp$hour == show[j]
      segments(0, 0, max(x[i]), sl[j] * max(x[i]), col = hcol[as.character(show[j])], lwd = 3)
    }
    stat_box(y, x, sprintf("within-site R² = %.3f",
                           within_r2(transform(cmp, pred = cmp[[p]]))))
    legend("bottomright", bty = "n", cex = 0.85, lty = c(3, rep(1, length(show))),
           lwd = c(1, rep(3, length(show))), col = c("grey40", hcol[as.character(show)]),
           legend = c("1:1", sprintf("%02d:00 slope = %.2f", show, sl)))
  }
  par(mar = c(6, 0.5, 6, 3))
  image(1, seq_along(hrs), t(seq_along(hrs)), col = hcol, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = seq_along(hrs), labels = sprintf("%02d:00", hrs), las = 1, cex.axis = 0.8)
  mtext("Hour", side = 3, line = 0.5, cex = 0.8)
  dev.off()
})

# ── Fig 9: daytime-integrated PAR ────────────────────────────────────────────
local({
  di <- lapply(names(CMP_MODELS), function(p) daytime_integrated(cmp, p))
  names(di) <- names(CMP_MODELS)
  ev_pch <- setNames(rep_len(c(0:14, 15:18), length(events)), events)
  dr     <- range(di[[1]]$depth)
  dcol   <- function(z) colorRampPalette(c("#FBD1B0", "#E8743B", "#2B0F06"))(100)[
                          1 + round(99 * (z - dr[1]) / diff(dr))]
  lim <- c(0, max(unlist(lapply(di, function(d) c(d$obs, d$pred)))))
  fig("fig9_DLI_M1_vs_M4r", w = 12, h = 6)
  layout(matrix(c(1, 2, 3), 1), widths = c(1, 1, 0.42))
  par(mar = c(4.5, 4.8, 3, 1))
  for (p in names(CMP_MODELS)) {
    d <- di[[p]]
    plot(d$pred, d$obs, pch = ev_pch[d$sitekey], col = dcol(d$depth), lwd = 1.3,
         xlim = lim, ylim = lim, main = CMP_MODELS[[p]], cex.main = 1,
         xlab = expression("Predicted daytime-integrated PAR (mol m"^-2 * " d"^-1 * ")"),
         ylab = expression("Measured daytime-integrated PAR (mol m"^-2 * " d"^-1 * ")"))
    abline(0, 1, lty = 2, col = "red", lwd = 1.5)
    abline(stats::lm(obs ~ pred, data = d), lwd = 2.5)
    fmt_stat_box(d$obs, d$pred)
  }
  par(mar = c(4, 0, 3, 0)); plot.new()
  legend("topleft", legend = events, pch = ev_pch, bty = "n", cex = 0.85,
         title = "Event", ncol = 2)
  zz <- seq(dr[1], dr[2], length.out = 6)
  legend("bottomleft", legend = sprintf("%.1f m", zz), pch = 15, col = dcol(zz),
         bty = "n", cex = 0.85, title = "Mean water depth")
  dev.off()
})

# ── S7: residual breakdown by event (pred - meas) ────────────────────────────
#  Points and boxes = primary model (M4r); red line = its loess with 95%
#  band; grey dashed = baseline (M1) loess. Open circles = binned means
#  ± SD (8 equal-count bins). Four events per page.
local({
  smooth_band <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 20 || length(unique(x[ok])) < 6) return(NULL)
    dd  <- data.frame(x = x[ok], y = y[ok])
    fit <- tryCatch(suppressWarnings(stats::loess(y ~ x, data = dd, span = 0.75)),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    g <- seq(min(dd$x), max(dd$x), length.out = 80)
    p <- tryCatch(suppressWarnings(stats::predict(fit, data.frame(x = g), se = TRUE)),
                  error = function(e) NULL)
    if (is.null(p) || all(!is.finite(p$fit))) return(NULL)
    data.frame(x = g, fit = p$fit, lo = p$fit - 1.96 * p$se.fit, hi = p$fit + 1.96 * p$se.fit)
  }
  bin_means <- function(x, y, n = 8) {
    br <- unique(stats::quantile(x, seq(0, 1, length.out = n + 1), na.rm = TRUE))
    if (length(br) < 3) return(NULL)
    b <- cut(x, br, include.lowest = TRUE)
    data.frame(x = tapply(x, b, mean), m = tapply(y, b, mean), s = tapply(y, b, stats::sd))
  }
  cmp$res_ref   <- cmp$pred_ref - cmp$light_mol
  cmp$res_base <- cmp$pred_base - cmp$light_mol
  ylim <- stats::quantile(c(cmp$res_ref, cmp$res_base), c(0.005, 0.995))
  pages <- split(events, ceiling(seq_along(events) / 4))
  for (pg in seq_along(pages)) {
    fig(sprintf("residuals_by_event/S7_residuals_by_event_p%d", pg), w = 13, h = 3.4 * length(pages[[pg]]) + 0.8)
    par(mfrow = c(length(pages[[pg]]), 3), mar = c(4.2, 4.5, 3.2, 1), oma = c(0, 0, 3, 0))
    for (k in pages[[pg]]) {
      d  <- cmp[cmp$sitekey == k, ]
      bx <- sprintf("bias: M4r %+.0f | M1 %+.0f  (n=%d)",
                    mean(d$res_ref), mean(d$res_base), nrow(d))
      boxplot(res_ref ~ hour, data = d, outline = FALSE, col = "#EAD9D5", border = "#7A5C55",
              ylim = ylim, xlab = "Hour of day", ylab = "Residual (pred - meas)")
      title(paste0(k, ": residual vs hour"), line = 1.4)
      med_c <- tapply(d$res_base, d$hour, stats::median)
      points(seq_along(med_c), med_c, pch = 4, col = "grey35", lwd = 2)
      abline(h = 0, lty = 2, col = "red"); mtext(bx, side = 3, line = 0.2, cex = 0.6)
      for (v in c("water_depth_m", "shadow_refr")) {
        plot(d[[v]], d$res_ref, pch = 16, cex = 0.4, col = adjustcolor("grey30", 0.35),
             ylim = ylim, ylab = "Residual (pred - meas)",
             xlab = if (v == "water_depth_m") "Water depth (m)" else
                      "Direct-light fraction (0 = shaded, 1 = lit)",
             main = "")
        title(paste0(k, ": residual vs ", if (v == "water_depth_m") "depth" else
                                            "direct-light fraction"), line = 1.4)
        abline(h = 0, lty = 2, col = "red")
        s <- smooth_band(d[[v]], d$res_ref)
        if (!is.null(s)) {
          polygon(c(s$x, rev(s$x)), c(s$lo, rev(s$hi)), col = adjustcolor("#C8372D", 0.2),
                  border = NA)
          lines(s$x, s$fit, col = "#C8372D", lwd = 2.5)
        }
        sc <- smooth_band(d[[v]], d$res_base)
        if (!is.null(sc)) lines(sc$x, sc$fit, col = "grey35", lwd = 2, lty = 2)
        b <- bin_means(d[[v]], d$res_ref)
        if (!is.null(b)) {
          segments(b$x, b$m - b$s, b$x, b$m + b$s, col = "grey20")
          points(b$x, b$m, pch = 21, bg = "white", cex = 1.2, lwd = 1.5)
        }
        mtext(bx, side = 3, line = 0.2, cex = 0.6)
      }
    }
    mtext(sprintf(paste("Residuals by event (%d of %d). Red = primary model (%s),",
                        "grey dashed / x = baseline M1d (constant Kd %.1f + direct shading)"),
                  pg, length(pages), FIG_MAIN_MODEL, KD_CONSTANT),
          outer = TRUE, line = 1, font = 2)
    dev.off()
  }
})

# ── Fig 10 + S8: event curves, baseline (M1) vs primary model (M4r) ──────────
# [NEW] Claude, 22 Oct 2026. Replaces the thesis Fig 8 grid in the paper.
#   One row per event, M1 left, M4r right, shared y axis within a row.
#   Grey = measured, orange = predicted: shaded band = range across loggers
#   (min-max) each hour, faint lines = individual loggers, bold line = mean.
#   x = local hour; a 24-h deployment means the morning hours come from the
#   day after the afternoon hours. Separate days are drawn as separate curves.
#   Fig 10: best, median and worst event by primary-model hourly R², chosen by
#     that rule in the code (or set FIG_EXAMPLE_EVENTS in config.R).
#   S8: all events, six per page, same layout.
#   The per-logger S1 time series above are still drawn as before.
local({
  ev_r2 <- sapply(events, function(k) {
    d <- cmp[cmp$sitekey == k, ]; metrics(d$light_mol, d$pred_ref)$r2 })
  ord <- names(sort(ev_r2, decreasing = TRUE))
  ex  <- if (exists("FIG_EXAMPLE_EVENTS") && length(FIG_EXAMPLE_EVENTS))
           FIG_EXAMPLE_EVENTS else ord[c(1, ceiling(length(ord) / 2), length(ord))]
  message("Fig 10 example events (best, median, worst by ", FIG_MAIN_MODEL, " R2): ",
          paste(ex, collapse = ", "))
  COL_MEAS <- "grey35"; COL_PRED <- "#E07B39"

  band_panel <- function(d, p, ylim, title, show_legend = FALSE) {
    d$pp  <- d[[p]]
    d$day <- format(d$datetime_pos, "%Y-%m-%d", tz = LOCAL_TZ)
    d$hr  <- as.numeric(format(d$datetime_pos, "%H", tz = LOCAL_TZ)) + 0.5
    plot(NA, xlim = c(6.5, 17.5), ylim = ylim, xaxt = "n",
         xlab = "Hour (local)", ylab = expression("PAR (" * mu * "mol m"^-2 * " s"^-1 * ")"))
    axis(1, at = seq(7, 17, 2))
    title(title, line = 0.6, cex.main = 0.95)
    for (dy in unique(d$day)) {
      s <- d[d$day == dy, ]
      for (v in c("light_mol", "pp")) {
        col <- if (v == "light_mol") COL_MEAS else COL_PRED
        a <- s %>% group_by(hr) %>%
          summarise(lo = min(.data[[v]]), hi = max(.data[[v]]), m = mean(.data[[v]]),
                    .groups = "drop") %>% arrange(hr)
        if (nrow(a) > 1)
          polygon(c(a$hr, rev(a$hr)), c(a$lo, rev(a$hi)), col = adjustcolor(col, 0.18), border = NA)
        for (id in unique(s$hoboID_env)) {
          l <- s[s$hoboID_env == id, ]; l <- l[order(l$hr), ]
          lines(l$hr, l[[v]], col = adjustcolor(col, 0.25), lwd = 0.6)
        }
        lines(a$hr, a$m, col = col, lwd = 3)
      }
    }
    m <- metrics(d$light_mol, d$pp)
    legend("topright", bty = "n", cex = 0.8,
           legend = c(sprintf("RMSE %.0f", m$rmse), sprintf("bias %+.0f", m$bias),
                      sprintf("R² %.2f", m$r2)))
    if (show_legend)
      legend("topleft", bty = "n", cex = 0.8, lwd = c(3, 3, 8, 8),
             col = c(COL_MEAS, COL_PRED, adjustcolor(COL_MEAS, 0.3), adjustcolor(COL_PRED, 0.3)),
             legend = c("measured, mean", "predicted, mean", "measured range", "predicted range"))
  }

  draw_rows <- function(keys, file, page_title) {
    fig(file, w = 10, h = 2.9 * length(keys) + 0.7)
    par(mfrow = c(length(keys), 2), mar = c(4, 4.5, 2.2, 1), oma = c(1.2, 2.5, 2.5, 0))
    for (k in keys) {
      d  <- cmp[cmp$sitekey == k, ]
      yl <- c(0, max(c(d$light_mol, d$pred_base, d$pred_ref)))
      band_panel(d, "pred_base", yl, "Baseline (M1d): constant Kd + direct shading",
                 show_legend = (k == keys[1]))
      mtext(k, side = 2, line = 5.2, font = 2, cex = 1.1)
      band_panel(d, "pred_ref", yl, sprintf("Primary (%s): reference Kd + shading", FIG_MAIN_MODEL))
    }
    mtext(page_title, outer = TRUE, line = 0.8, font = 2)
    dev.off()
  }

  draw_rows(ex, "fig10_event_examples",
            sprintf("Best, median and worst events by %s hourly R²", FIG_MAIN_MODEL))
  pages <- split(ord, ceiling(seq_along(ord) / 6))
  for (pg in seq_along(pages))
    draw_rows(pages[[pg]], sprintf("event_curves/S8_event_curves_p%d", pg),
              sprintf("All events, ordered by %s hourly R² (page %d of %d)",
                      FIG_MAIN_MODEL, pg, length(pages)))
})

write_manifest("04_figures")
message("04 done. Figures in ", out_path("figures"))
