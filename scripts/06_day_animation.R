# ══════════════════════════════════════════════════════════════════════════════
#  06_day_animation.R: ray-trace one site across a whole day
#
#  Usage:   Rscript scripts/06_day_animation.R            (site = FIG_SITE)
#           Rscript scripts/06_day_animation.R CB18 10    (site, minutes per step)
#
#  Reads : output/dems/<SITE>.rds, output/geom.rds, output/cache.rds
#  Writes: output/figures/day_<SITE>/frame_###.png   one frame per time step
#          output/figures/day_<SITE>.gif              if gifski or magick is installed
#          output/figures/day_<SITE>_lit_fraction.png fraction of the day each
#                                                    pixel is in direct sun,
#                                                    refracted vs unrefracted
#          output/figures/day_<SITE>_stack.rds        the shadow stacks
#
#  Uses the same shadow function as stage 02 (make_shadow), so it inherits the
#  orientation and sun-direction corrections. Purely a visual check: nothing
#  here feeds the models.
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_raytrace.R"))
suppressPackageStartupMessages({ library(dplyr) })   # load_hobo_locations() uses %>%

args <- commandArgs(trailingOnly = TRUE)
SITE <- if (length(args) >= 1) args[1] else FIG_SITE
STEP <- if (length(args) >= 2) as.numeric(args[2]) else 15   # minutes

stage_banner(sprintf("06  DAY ANIMATION  (%s, every %g min)", SITE, STEP))

# Same orientation / sun-direction set-up as stage 02. The messages printed
# here must match stage 02's (e.g. "using 'transpose'", "s = +1, offset = 0°").
init_geometry_checks()

dem   <- load_dem(in_path("dems", paste0(SITE, ".rds")))
hobo  <- load_hobo_locations(readRDS(in_path("geom.rds")), SITE)
cache <- readRDS(in_path("cache.rds"))
cen   <- get_site_centroid(dem$rast)

# the day with the most logger-hours at this site
ts  <- cache$datetime_pos[cache$sitekey == SITE]
day <- as.Date(names(sort(table(as.Date(ts, tz = LOCAL_TZ)), decreasing = TRUE))[1])

times <- seq(as.POSIXct(paste(day, "05:00"), tz = LOCAL_TZ),
             as.POSIXct(paste(day, "19:00"), tz = LOCAL_TZ), by = STEP * 60)

# sun_at_midhour() evaluates the sun 30 min after the time it is given,
# so pass t - 30 min to get the sun at t itself.
sun_at <- function(t) sun_at_midhour(t - 1800, cen$lat, cen$lon)

alts  <- vapply(times, function(t) sun_at(t)$altitude * 180 / pi, numeric(1))
times <- times[alts > SUN_ALT_MIN_DEG]
message(sprintf("  %s on %s: %d time steps from %s to %s", SITE, day, length(times),
                format(min(times), "%H:%M"), format(max(times), "%H:%M")))

frame_dir <- out_path("figures", paste0("day_", SITE))
dir.create(frame_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(frame_dir, "frame_*.png"))

# mask off-reef pixels on the matrix, exactly as stage 02 does, then convert
shade_rast <- function(sun, refract) {
  m <- make_shadow(dem, sun, refract = refract)
  m[!dem$reef_mask] <- NA
  matrix_to_rast(m, dem$rast)
}
pal <- gray.colors(100, start = 0.1, end = 0.98)
fig <- function(name, w = 8, h = 6)   # same helper as 04_figures.R
  png(out_path("figures", paste0(name, ".png")), width = w, height = h, units = "in", res = 200)

stack_refr <- list(); stack_air <- list()
for (i in seq_along(times)) {
  sun <- sun_at(times[i])
  r_air  <- shade_rast(sun, refract = FALSE)
  r_refr <- shade_rast(sun, refract = TRUE)
  stack_air[[i]] <- r_air; stack_refr[[i]] <- r_refr

  png(file.path(frame_dir, sprintf("frame_%03d.png", i)), width = 1600, height = 820, res = 110)
  par(mfrow = c(1, 2), mar = c(0.5, 0.5, 3, 0.5), oma = c(0, 0, 2.5, 0))
  for (k in 1:2) {
    r <- terra::trim(if (k == 1) r_air else r_refr)
    a <- raytrace_angles(sun, refract = (k == 2))
    terra::plot(r, col = pal, range = c(0, 1), legend = FALSE, axes = FALSE,
                main = sprintf("%s  sun %.0f°", if (k == 1) "Unrefracted (in air)" else "Refracted (in water)", a$alt_deg))
    terra::plot(terra::vect(hobo), add = TRUE, pch = 21, col = "red", cex = 0.9)
  }
  mtext(sprintf("%s  %s %s   azimuth %.0f°", SITE, day, format(times[i], "%H:%M"),
                (sun$azimuth * 180 / pi) %% 360), outer = TRUE, line = 0.5, cex = 1.3)
  dev.off()
  message(sprintf("  frame %3d/%d  %s", i, length(times), format(times[i], "%H:%M")))
}

# ── Fraction of the day each pixel is in direct sun ─────────────────────────
lit_air  <- terra::mean(terra::rast(stack_air))
lit_refr <- terra::mean(terra::rast(stack_refr))
names(lit_air) <- "unrefracted"; names(lit_refr) <- "refracted"

fig(paste0("day_", SITE, "_lit_fraction"), w = 13, h = 6.5)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 3, 4), oma = c(0, 0, 2, 0))
for (r in list(lit_air, lit_refr)) {
  terra::plot(terra::trim(r), col = hcl.colors(100, "viridis"), range = c(0, 1),
              axes = FALSE, main = names(r))
  terra::plot(terra::vect(hobo), add = TRUE, pch = 21, col = "red", cex = 0.9)
}
mtext(sprintf("%s %s: fraction of daytime (sun > %g°) in direct sun", SITE, day, SUN_ALT_MIN_DEG),
      outer = TRUE, line = 0.3)
dev.off()

v_air <- terra::values(lit_air)[, 1]; v_refr <- terra::values(lit_refr)[, 1]
message(sprintf("  reef pixels fully lit all day: unrefracted %.0f%%, refracted %.0f%%",
                100 * mean(v_air >= 0.999, na.rm = TRUE), 100 * mean(v_refr >= 0.999, na.rm = TRUE)))
message(sprintf("  SD of daily lit fraction:     unrefracted %.3f, refracted %.3f",
                sd(v_air, na.rm = TRUE), sd(v_refr, na.rm = TRUE)))

saveRDS(list(times = times, refracted = terra::wrap(terra::rast(stack_refr)),
                    unrefracted = terra::wrap(terra::rast(stack_air))),
               out_path("figures", paste0("day_", SITE, "_stack.rds")))

# ── GIF ─────────────────────────────────────────────────────────────────────
frames <- sort(list.files(frame_dir, "^frame_.*\\.png$", full.names = TRUE))
gif    <- out_path("figures", paste0("day_", SITE, ".gif"))
if (requireNamespace("gifski", quietly = TRUE)) {
  gifski::gifski(frames, gif, width = 1600, height = 820, delay = 0.25, progress = FALSE)
  message("  GIF: ", gif)
} else if (requireNamespace("magick", quietly = TRUE)) {
  magick::image_write(magick::image_animate(magick::image_read(frames), fps = 4), gif)
  message("  GIF: ", gif)
} else {
  message("  No gifski/magick: frames are in ", frame_dir,
          " (install.packages('gifski') for a GIF)")
}
message("06 done.")
