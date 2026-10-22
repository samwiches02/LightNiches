# ══════════════════════════════════════════════════════════════════════════════
#  08_logger_check.R — is this logger's shading being read from the right place?
#
#  Re-runs the ray-trace for ONE event at chosen hours, re-extracts the shading
#  at one logger, and compares it with the value stage 02 stored in cache.rds.
#  Run on demand; not part of run_all.R. Needs 02 (and ideally 03) outputs.
#
#  From a terminal in the project folder:
#     Rscript scripts/08_logger_check.R HS19 H14
#         -> 08:00, 12:00, 16:00 (FIG_HOURS) on the day with the most data
#     Rscript scripts/08_logger_check.R HS19 H14 2019-11-10
#     Rscript scripts/08_logger_check.R HS19 H14 2019-11-10 9,13
#     Rscript scripts/08_logger_check.R L119 H05,H22,H31 2019-11-04 10
#         -> one figure per logger
#  From RStudio:
#     source("scripts/08_logger_check.R")
#     logger_check("HS19", "H14", day = "2019-11-10", hours = c(8, 12, 16))
#
#  Writes: output/figures/logger_check/<event>_<logger>_<day>_h<hours>.png
#
#  One column per hour:
#    row 1  whole-reef direct-light map; red = all loggers, yellow = this one,
#           blue square = the zoom window
#    row 2  measured vs predicted PAR for every logger that hour (this one boxed)
#    row 3  zoom on the logger with the extraction buffer; title gives the value
#           re-extracted now and the value stored by stage 02
#  Row 4 (once): elevation zoom, sky-view (ambient) zoom, and DEM elevation vs
#    dpHB for every logger on the event (stage 05's registration test). If the
#    positions are registered, z tracks dpHB; this logger's distance from the
#    line says whether its position or its depth is off.
#
#  How to read it:
#    "now" = "stage 02" (MATCH)  -> the cached shading is what this map gives.
#    MISMATCH                    -> the cache came from a different DEM or
#                                   positions (e.g. re-registered DEM, stale
#                                   output/raytrace/<event>.rds). Re-run 02.
#    MATCH but the zoom shows the logger beside an obvious shadow edge, or the
#    bars show a lit logger measured dark (or the reverse) -> position, DEM
#    detail or an overhang the top-down DEM cannot represent.
#
#  PROVENANCE: [NEW] Claude, 22 Oct 2026. Layout follows Sam's
#    plot_shadow_3hr_with_zoom(): three hours, whole-reef map, bar chart, zoom.
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_raytrace.R"))
suppressPackageStartupMessages({ library(dplyr) })

LOGGER_CHECK_ZOOM_M <- 0.25    # half-width of the zoom window (m)
LOGGER_CHECK_TOL    <- 0.01    # |now − stage 02| above this is flagged MISMATCH

logger_check <- function(key, logger, day = NULL, hours = FIG_HOURS,
                         model = FIG_MAIN_MODEL, buffer_m = BUFFER_M,
                         zoom_m = LOGGER_CHECK_ZOOM_M, show = interactive()) {
  if (!exists(".SHADE_FIX", envir = globalenv())) init_geometry_checks()

  # ── inputs ──
  dem_file <- in_path("dems", paste0(key, ".rds"))
  if (!file.exists(dem_file)) stop("No prepared DEM for ", key, " (", dem_file,
                                   "). Run stage 02 first.", call. = FALSE)
  dem   <- load_dem(dem_file)
  hobo  <- load_hobo_locations(readRDS(in_path("geom.rds")), key)
  if (!logger %in% hobo$hoboID)
    stop(logger, " has no position on ", key, ". Loggers with positions: ",
         paste(sort(hobo$hoboID), collapse = " "), call. = FALSE)
  cache <- readRDS(in_path("cache.rds")) %>% filter(sitekey == key)
  preds <- if (file.exists(in_path("predictions.rds")))
    readRDS(in_path("predictions.rds")) %>% filter(sitekey == key, model == !!model)
  else NULL

  # ── which day and hours ──
  mine <- cache$datetime_pos[cache$hoboID_env == logger]
  if (length(mine) == 0)
    message("  Note: ", logger, " has a position but no rows in cache.rds ",
            "(no light data, or it was dropped in 01/02). Shading is still drawn.")
  pool <- if (length(mine)) mine else cache$datetime_pos
  if (is.null(day))
    day <- names(sort(table(as.Date(pool, tz = LOCAL_TZ)), decreasing = TRUE))[1]
  day   <- as.Date(day)
  times <- as.POSIXct(sprintf("%s %02d:00:00", day, as.integer(hours)), tz = LOCAL_TZ)

  # ── geometry that does not depend on the sun ──
  cen    <- get_site_centroid(dem$rast)
  me     <- hobo[hobo$hoboID == logger, ]
  xy     <- sf::st_coordinates(me)[1, ]
  win    <- terra::ext(xy[1] - zoom_m, xy[1] + zoom_m, xy[2] - zoom_m, xy[2] + zoom_m)
  circle <- terra::vect(sf::st_buffer(me, buffer_m))
  others <- terra::vect(hobo)
  amb    <- dem$ambient; amb[!dem$reef_mask] <- NA
  elev   <- terra::mask(dem$rast, matrix_to_rast(dem$reef_mask * 1, dem$rast),
                        maskvalues = 0)
  pal    <- gray.colors(100, start = 0.05, end = 0.98)
  cached <- function(t, col) {
    v <- cache[[col]][cache$hoboID_env == logger & cache$datetime_pos == t]
    if (length(v)) v[1] else NA_real_
  }
  verdict <- function(now, then) {
    if (!is.finite(then)) return("no stage-02 value")
    if (abs(now - then) <= LOGGER_CHECK_TOL) "MATCH" else "MISMATCH"
  }

  # all loggers: DEM elevation vs dpHB (stage 05's registration test)
  z_all <- terra::extract(dem$rast, terra::vect(sf::st_buffer(hobo, buffer_m)),
                          fun = mean, na.rm = TRUE)[, 2]
  dp    <- cache %>% distinct(hoboID_env, dpHB)
  reg   <- data.frame(hoboID = hobo$hoboID, z = z_all) %>%
    inner_join(dp, by = c("hoboID" = "hoboID_env"))
  amb_now <- extract_at_loggers(amb, dem$rast, me, buffer_m)

  # ── figure ──
  nh  <- length(times)
  dir <- out_path("figures", "logger_check", "x"); dir <- dirname(dir)
  png_file <- file.path(dir, sprintf("%s_%s_%s_h%s.png", key, logger, day,
                                      paste(sprintf("%02d", as.integer(hours)), collapse = "-")))
  png(png_file, width = 4.6 * max(nh, 3), height = 16, units = "in", res = 150)
  nc <- max(nh, 3)
  lay <- rbind(seq_len(nc), nc + seq_len(nc), 2 * nc + seq_len(nc), 3 * nc + seq_len(nc))
  layout(lay, heights = c(1.15, 0.7, 1, 1))
  par(oma = c(0, 0, 3.5, 0))
  maps <- vector("list", nh); vals <- numeric(nh)

  # row 1: whole-reef map
  for (i in seq_len(nh)) {
    sun <- sun_at_midhour(times[i], cen$lat, cen$lon)
    a   <- raytrace_angles(sun, refract = TRUE)
    maps[[i]] <- matrix_to_rast(make_shadow(dem, sun, refract = TRUE), dem$rast)
    vals[i]   <- terra::extract(maps[[i]], circle, fun = mean, na.rm = TRUE)[1, 2]  # as in 02
    par(mar = c(0.5, 0.5, 3, 0.5))
    if (sun$altitude * 180 / pi <= SUN_ALT_MIN_DEG) {
      plot.new(); title(sprintf("%s  sun %.0f° in air: below SUN_ALT_MIN_DEG, not modelled",
                                format(times[i], "%H:%M"), sun$altitude * 180 / pi))
      next
    }
    terra::plot(terra::trim(maps[[i]]), col = pal, range = c(0, 1), legend = FALSE,
                axes = FALSE, mar = NA)
    terra::plot(others, add = TRUE, pch = 21, col = "red", cex = 1.1, lwd = 1.3)
    terra::plot(terra::as.polygons(win), add = TRUE, border = "dodgerblue", lwd = 2)
    points(xy[1], xy[2], pch = 23, bg = "yellow", col = "black", cex = 1.8)
    title(sprintf("%s  sun alt %.0f° (in water %.0f°), bearing %.0f°",
                  format(times[i], "%H:%M"), sun$altitude * 180 / pi, a$alt_deg, a$bearing),
          line = 0.8)
  }
  for (i in seq_len(nc - nh)) plot.new()

  # row 2: bars, every logger at that hour
  for (i in seq_len(nh)) {
    par(mar = c(4, 4, 1.5, 0.5))
    p <- if (!is.null(preds)) preds[preds$datetime_pos == times[i], ] else NULL
    if (is.null(p) || nrow(p) == 0) {
      plot.new(); title(sprintf("no %s predictions at %s", model, format(times[i], "%H:%M")),
                        line = 0); next
    }
    p   <- p[order(p$hoboID_env), ]
    yl  <- c(0, max(c(p$light_mol, p$pred), na.rm = TRUE) * 1.08)
    bp  <- barplot(rbind(p$pred, p$light_mol), beside = TRUE, col = c("firebrick", "steelblue"),
                   border = NA, names.arg = p$hoboID_env, las = 2, cex.names = 0.6, ylim = yl,
                   ylab = expression("PAR (" * mu * "mol m"^-2 * " s"^-1 * ")"))
    j <- match(logger, p$hoboID_env)
    if (!is.na(j)) rect(bp[1, j] - 0.6, 0, bp[2, j] + 0.6, yl[2], border = "dodgerblue", lwd = 2)
    if (i == 1) legend("topleft", c(paste("Predicted", model), "Measured"), bty = "n",
                       fill = c("firebrick", "steelblue"), border = NA, cex = 0.85)
  }
  for (i in seq_len(nc - nh)) plot.new()

  # row 3: zoom on the logger, direct light
  for (i in seq_len(nh)) {
    par(mar = c(3, 3, 3, 0.5))
    then <- cached(times[i], "shadow_refr")
    terra::plot(terra::crop(maps[[i]], win), col = pal, range = c(0, 1), legend = FALSE,
                mar = NA, cex.axis = 0.6)
    terra::plot(circle, add = TRUE, border = "dodgerblue", col = adjustcolor("dodgerblue", 0.2),
                lwd = 2)
    points(xy[1], xy[2], pch = 23, bg = "yellow", cex = 1.6)
    title(sprintf("direct light: now %.2f | stage 02 %s  [%s]", vals[i],
                  if (is.finite(then)) sprintf("%.2f", then) else "-",
                  verdict(vals[i], then)),
          line = 0.8, col.main = if (identical(verdict(vals[i], then), "MISMATCH")) "red" else "black")
  }
  for (i in seq_len(nc - nh)) plot.new()

  # row 4: elevation zoom, sky-view zoom, registration scatter
  par(mar = c(3, 3, 3, 0.5))
  terra::plot(terra::crop(elev, win), col = hcl.colors(100, "Terrain"), mar = NA,
              cex.axis = 0.6, plg = list(cex = 0.7))
  terra::plot(circle, add = TRUE, border = "black", lwd = 2)
  points(xy[1], xy[2], pch = 23, bg = "yellow", cex = 1.6)
  z_me <- reg$z[reg$hoboID == logger]; dp_me <- reg$dpHB[reg$hoboID == logger]
  title(sprintf("DEM elevation: %.2f m at logger | dpHB %s", z_me[1],
                if (length(dp_me)) sprintf("%.2f m", dp_me) else "-"), line = 0.8)

  terra::plot(terra::crop(matrix_to_rast(amb, dem$rast), win), col = pal, range = c(0, 1),
              legend = FALSE, mar = NA, cex.axis = 0.6)
  terra::plot(circle, add = TRUE, border = "dodgerblue", lwd = 2)
  points(xy[1], xy[2], pch = 23, bg = "yellow", cex = 1.6)
  amb_then <- if (length(mine)) cached(mine[1], "ambient") else NA_real_
  title(sprintf("sky view: now %.2f | stage 02 %s  [%s]", amb_now,
                if (is.finite(amb_then)) sprintf("%.2f", amb_then) else "-",
                verdict(amb_now, amb_then)), line = 0.8)

  par(mar = c(4.5, 4.5, 3, 1))
  if (nrow(reg) >= 3) {
    fit <- stats::lm(z ~ dpHB, data = reg)
    plot(reg$dpHB, reg$z, pch = 21, bg = "grey80", xlab = "dpHB (m)",
         ylab = "DEM elevation at logger (m)")
    abline(fit, lty = 2)
    k <- reg$hoboID == logger
    points(reg$dpHB[k], reg$z[k], pch = 23, bg = "yellow", cex = 2)
    text(reg$dpHB, reg$z, reg$hoboID, cex = 0.5, pos = 4, col = "grey40")
    res_me <- if (any(k)) stats::residuals(fit)[k] else NA
    title(sprintf("registration: slope %.2f, r %.2f | %s residual %+.2f m",
                  stats::coef(fit)[2], stats::cor(reg$dpHB, reg$z), logger, res_me),
          line = 0.8)
  } else { plot.new(); title("too few loggers with depth for the registration plot") }
  for (i in seq_len(nc - 3)) plot.new()

  mtext(sprintf("%s  |  logger %s  |  %s  |  buffer %.0f cm  |  model %s",
                key, logger, format(day, "%d %b %Y"), buffer_m * 100, model),
        outer = TRUE, line = 1.2, font = 2, cex = 1.2)
  dev.off()
  message("  Wrote ", png_file)

  out <- data.frame(sitekey = key, logger = logger, time = format(times, "%Y-%m-%d %H:%M"),
                    shadow_now = vals,
                    shadow_stage02 = vapply(times, function(t) cached(t, "shadow_refr"), 0),
                    ambient_now = amb_now, ambient_stage02 = amb_then)
  out$verdict <- mapply(verdict, out$shadow_now, out$shadow_stage02)
  print(out, row.names = FALSE)
  if (show && interactive()) utils::browseURL(png_file)
  invisible(out)
}

# ── command line ──
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2)
    stop("Usage: Rscript scripts/08_logger_check.R <event> <logger[,logger...]> ",
         "[YYYY-MM-DD] [hours, e.g. 8,12,16]", call. = FALSE)
  day   <- if (length(args) >= 3) args[3] else NULL
  hours <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else FIG_HOURS
  for (lg in strsplit(args[2], ",")[[1]]) logger_check(args[1], lg, day, hours)
}
