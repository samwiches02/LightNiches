# ══════════════════════════════════════════════════════════════════════════════
#  05_diagnostics.R — are logger positions registered correctly to the DEM?
#
#  PROVENANCE: [NEW] — written by Claude. Not part of run_all.R; run by hand
#  after stage 02 (it does not need stage 03).
#
#  IDEA
#  Most of the measured logger-to-logger spread is a constant per-logger offset
#  (some loggers read high all day). Sky view (the ambient layer) is also
#  constant through the day, so if positions are right, measured offsets should
#  correlate with sky view at the logger locations. If the positions are
#  rotated, mirrored or shifted relative to the DEM, the correlation shows up
#  under the matching transformation instead.
#
#  No ray-tracing: it re-extracts the saved ambient layer from output/dems/.
#
#  Reads : output/cache.rds, output/geom.rds, output/dems/<SITEyy>.rds
#  Writes: output/diagnostics/registration_<SITEyy>.csv
#          output/diagnostics/shift_<SITEyy>.csv
#          output/diagnostics/buffer_<SITEyy>.csv
#
#  Run:  Rscript scripts/05_diagnostics.R      (or open and Source in RStudio)
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_raytrace.R"))
suppressPackageStartupMessages({ library(dplyr) })

stage_banner("05  REGISTRATION DIAGNOSTICS")
init_geometry_checks()
set.seed(SEED)

cache <- readRDS(out_path("cache.rds"))
geom  <- readRDS(out_path("geom.rds"))

# ── Measured per-logger offset: mean log brightness relative to the site-hour mean
offsets <- cache %>%
  filter(is.finite(light_mol), light_mol > 0) %>%
  group_by(sitekey, datetime_pos) %>% filter(dplyr::n() >= 3) %>%
  mutate(lo = log(light_mol / mean(light_mol))) %>%
  group_by(sitekey, hoboID_env) %>%
  summarise(offset = mean(lo), depth = mean(water_depth_m), n_hours = dplyr::n(),
            .groups = "drop") %>%
  filter(n_hours >= 3)

# ── Helpers ──────────────────────────────────────────────────────────────────
# Rotate (deg, anticlockwise), optionally mirror, then shift points about a centre.
transform_pts <- function(xy, centre, rot = 0, mirror = FALSE, dx = 0, dy = 0) {
  x <- xy[, 1] - centre[1]; y <- xy[, 2] - centre[2]
  if (mirror) x <- -x                                   # mirror east-west
  a  <- rot * pi / 180
  xr <- x * cos(a) - y * sin(a); yr <- x * sin(a) + y * cos(a)
  cbind(xr + centre[1] + dx, yr + centre[2] + dy)
}

ambient_at <- function(amb_r, xy, buffer_m) {
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"),
                      crs = UTM_CRS)
  terra::extract(amb_r, terra::vect(sf::st_buffer(pts, buffer_m)),
                 fun = mean, na.rm = TRUE)[, 2]
}

scor <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman",
                                                  use = "complete.obs"))

# ── Per site ─────────────────────────────────────────────────────────────────
summary_rows <- list()
for (key in SITES) {
  f <- out_path("dems", paste0(key, ".rds"))
  if (!file.exists(f)) { message(key, ": no saved DEM, skipping"); next }
  message("\n── ", key, " ──")
  dem <- load_dem(f)
  amb <- dem$ambient; amb[!dem$reef_mask] <- NA
  amb_r <- matrix_to_rast(amb, dem$rast)
  centre <- colMeans(terra::xyFromCell(amb_r, which(!is.na(terra::values(amb_r)))))

  hobo <- load_hobo_locations(geom, key)
  off  <- offsets %>% filter(sitekey == key)
  hobo <- hobo[match(off$hoboID_env, hobo$hoboID), ]
  ok   <- !is.na(hobo$hoboID)
  off  <- off[ok, ]; hobo <- hobo[ok, ]
  xy   <- sf::st_coordinates(hobo)
  message("  ", nrow(off), " loggers with offsets and positions")

  # 1) the 8 rotations x mirror, no shift
  reg <- expand.grid(rot = c(0, 90, 180, 270), mirror = c(FALSE, TRUE))
  reg$cor_ambient <- mapply(function(r, m)
    scor(off$offset, ambient_at(amb_r, transform_pts(xy, centre, r, m), BUFFER_M)),
    reg$rot, reg$mirror)
  # null: correlation with logger labels shuffled (positions random w.r.t. offsets)
  amb0 <- ambient_at(amb_r, xy, BUFFER_M)
  null <- replicate(1000, scor(sample(off$offset), amb0))
  reg$null_95 <- unname(stats::quantile(abs(null), 0.95, na.rm = TRUE))
  reg$label <- ifelse(reg$rot == 0 & !reg$mirror, "AS GIVEN", "")
  print(reg, digits = 2, row.names = FALSE)
  write.csv(reg, out_path("diagnostics", paste0("registration_", key, ".csv")), row.names = FALSE)

  # 2) shifts of the positions as given, -1 m to +1 m
  sh <- expand.grid(dx = seq(-1, 1, 0.1), dy = seq(-1, 1, 0.1))
  sh$cor_ambient <- mapply(function(a, b)
    scor(off$offset, ambient_at(amb_r, transform_pts(xy, centre, 0, FALSE, a, b), BUFFER_M)),
    sh$dx, sh$dy)
  best <- sh[which.max(sh$cor_ambient), ]
  write.csv(sh, out_path("diagnostics", paste0("shift_", key, ".csv")), row.names = FALSE)

  # 3) buffer sweep at the positions as given
  bs <- data.frame(buffer_cm = c(1, 2, 4, 6, 10, 20, 50))
  bs$cor_ambient <- vapply(bs$buffer_cm, function(b)
    scor(off$offset, ambient_at(amb_r, xy, b / 100)), numeric(1))
  write.csv(bs, out_path("diagnostics", paste0("buffer_", key, ".csv")), row.names = FALSE)

  summary_rows[[key]] <- data.frame(
    sitekey = key, n = nrow(off),
    cor_depth = scor(off$offset, -off$depth),               # + means deeper = darker
    cor_ambient_as_given = reg$cor_ambient[reg$label == "AS GIVEN"],
    best_transform = with(reg[which.max(reg$cor_ambient), ],
                          sprintf("rot %d%s", rot, if (mirror) " + mirror" else "")),
    best_transform_cor = max(reg$cor_ambient, na.rm = TRUE),
    best_shift = sprintf("dx %.1f, dy %.1f m", best$dx, best$dy),
    best_shift_cor = best$cor_ambient,
    null_95 = reg$null_95[1],
    best_buffer_cm = bs$buffer_cm[which.max(bs$cor_ambient)])
}

res <- bind_rows(summary_rows)

# ══════════════════════════════════════════════════════════════════════════════
#  PART B — checks that do NOT use rayshader
#  [Claude] dpHB was read from the DEM at the logger positions in geometry.R.
#  If our positions line up with our DEM, DEM elevation at each logger should
#  match dpHB closely AS GIVEN. Relief (logger height above the mean height
#  within 50 cm) is a rayshader-free proxy for openness.
# ══════════════════════════════════════════════════════════════════════════════
message("\n\nPART B: DEM elevation vs dpHB, and relief vs measured offsets (no rayshader)")
partB <- list()
for (key in SITES) {
  f <- out_path("dems", paste0(key, ".rds")); if (!file.exists(f)) next
  dem <- load_dem(f)
  zr  <- dem$rast
  inside <- terra::values(matrix_to_rast(dem$reef_mask * 1, dem$rast))[, 1] > 0.5
  v   <- terra::values(zr)[, 1]; v[!inside] <- NA            # keep the reef disc only
  terra::values(zr) <- v
  centre <- colMeans(terra::xyFromCell(zr, which(!is.na(v))))
  hobo <- load_hobo_locations(geom, key)
  lg <- cache %>% filter(sitekey == key) %>%
    distinct(hoboID_env, dpHB) %>%
    left_join(offsets %>% filter(sitekey == key) %>% select(hoboID_env, offset), by = "hoboID_env")
  hobo <- hobo[match(lg$hoboID_env, hobo$hoboID), ]
  ok <- !is.na(hobo$hoboID); lg <- lg[ok, ]; hobo <- hobo[ok, ]
  xy <- sf::st_coordinates(hobo)
  tr <- expand.grid(rot = c(0, 90, 180, 270), mirror = c(FALSE, TRUE))
  out <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    p <- transform_pts(xy, centre, tr$rot[i], tr$mirror[i])
    z_pt  <- ambient_at(zr, p, 0.02)
    z_nb  <- ambient_at(zr, p, 0.50)
    data.frame(rot = tr$rot[i], mirror = tr$mirror[i],
               cor_z_dpHB = suppressWarnings(stats::cor(z_pt, lg$dpHB, use = "complete.obs")),
               cor_relief_offset = scor(z_pt - z_nb, lg$offset),
               label = if (tr$rot[i] == 0 && !tr$mirror[i]) "AS GIVEN" else "")
  }))
  message("\n── ", key, " ──")
  print(out, digits = 2, row.names = FALSE)
  partB[[key]] <- cbind(sitekey = key, out)
}
write.csv(bind_rows(partB), out_path("diagnostics", "partB_dem_checks.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  PART C — does rayshader put its OUTPUT where the terrain is?
#  [Claude] Off-centre pillar on a NON-SQUARE surface (a centred pillar cannot
#  reveal a reflected output). Sun due east via the calibration. The shadow
#  must touch the pillar on its WEST side, and the lowest sky view must be at
#  the pillar's foot, not at its mirror position.
# ══════════════════════════════════════════════════════════════════════════════
message("\n\nPART C: rayshader output alignment (off-centre pillar, non-square grid)")
local({
  r <- terra::rast(nrows = 40, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 40,
                   vals = 0, crs = UTM_CRS)
  r[6:10, 8:12] <- 5                                     # pillar near the NW corner
  m <- dem_to_matrix(r)
  pil <- colMeans(terra::xyFromCell(r, which(terra::values(r)[, 1] > 0)))
  s <- shade_fix(rayshader::ray_shade(m, sunaltitude = 30,
                            sunangle = sunangle_for_bearing(90, get(".SUN_CAL", envir = globalenv())),
                            zscale = 1, lambert = FALSE, progbar = FALSE))
  message(sprintf("  ray_shade output (after correction) dims %s vs input %s",
                  paste(dim(s), collapse = "x"), paste(dim(m), collapse = "x")))
  sr <- matrix_to_rast(s, r)
  g  <- which(terra::values(r)[, 1] == 0 & terra::values(sr)[, 1] < 0.5)
  sh <- colMeans(terra::xyFromCell(r, g))
  message(sprintf("  pillar centre (%.1f, %.1f); shadow centre (%.1f, %.1f); shadow is %.1f m %s of pillar, %.1f m %s",
                  pil[1], pil[2], sh[1], sh[2], abs(sh[1] - pil[1]),
                  if (sh[1] < pil[1]) "WEST" else "EAST", abs(sh[2] - pil[2]),
                  if (sh[2] < pil[2]) "south" else "north"))
  okC1 <- sh[1] < pil[1] && abs(sh[2] - pil[2]) < 3 && sqrt(sum((sh - pil)^2)) < 15
  message("  -> ", if (okC1) "OK: shadow is beside the pillar, on its west side"
                   else "PROBLEM: shadow is not beside the pillar on its west side")

  a  <- amb_fix(rayshader::ambient_shade(m, sunbreaks = 24, maxsearch = 20, zscale = 1,
                                         multicore = FALSE, progbar = FALSE))
  ar <- matrix_to_rast(a, r)
  gv <- terra::values(r)[, 1] == 0
  xyg <- terra::xyFromCell(r, which(gv))
  d  <- sqrt((xyg[, 1] - pil[1])^2 + (xyg[, 2] - pil[2])^2)
  rho <- scor(d, terra::values(ar)[gv, 1])
  message(sprintf("  sky view vs distance from pillar: Spearman %.2f -> %s", rho,
                  if (rho > 0.3) "OK: most enclosed ground is at the pillar's foot"
                  else "PROBLEM: sky-view map is not aligned with the terrain"))
})
message("\nSUMMARY  (correlations are Spearman, measured per-logger offset vs sky view)")
print(res, digits = 2, row.names = FALSE)
write.csv(res, out_path("diagnostics", "summary.csv"), row.names = FALSE)

message("\nHow to read it:",
        "\n  cor_ambient_as_given clearly positive and above null_95 -> positions OK;",
        "\n     offsets partly reflect real sky view.",
        "\n  A rotation/mirror much higher than AS GIVEN (and above null_95) -> positions",
        "\n     are rotated/mirrored relative to the DEM.",
        "\n  Best shift far from 0 with a much higher correlation -> positions are offset.",
        "\n  Everything inside +/- null_95 -> no detectable link between offsets and",
        "\n     geometry; sensor calibration is the likelier cause.")
