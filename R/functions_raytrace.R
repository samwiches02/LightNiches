# ══════════════════════════════════════════════════════════════════════════════
#  functions_raytrace.R — DEM preparation, sun position, shading, extraction
#  Used by scripts/02_raytrace.R and scripts/04_figures.R.
#
#  PROVENANCE TAGS (every function has one):
#    [YOURS]          your code from master_script.R, logic unchanged
#    [YOURS, EDITED]  your code, with the listed changes
#    [NEW]            written by Claude
#  Comments with no "[Claude]" marker are yours, copied from master_script.R.
# ══════════════════════════════════════════════════════════════════════════════


# ── Raster <-> rayshader matrix, with the orientation MEASURED, not assumed ───
# [NEW]
#   [Claude] HISTORY: master_script.R (and the first version of this file)
#   converted rayshader matrices back to rasters with rast(t(m)). The CB test
#   run showed this TRANSPOSES the map: logger offsets matched sky view only
#   after a transpose of the positions (05_diagnostics.R), and the sun-angle
#   calibration came out as rayshader's documented convention seen through a
#   transpose. The symmetric pillar test could not detect it.
#
#   Fix: convert a small NON-SQUARE raster with all-different values through
#   rayshader::raster_to_matrix, find which of the 8 rotations/flips gives the
#   original back exactly, and use that inverse everywhere. Any version of
#   rayshader is handled, and a mismatch stops the run.
dem_to_matrix <- function(r) rayshader::raster_to_matrix(raster::raster(r))

detect_matrix_orientation <- function() {
  r <- terra::rast(nrows = 7, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 7,
                   vals = 1:35, crs = UTM_CRS)
  A <- terra::as.matrix(r, wide = TRUE)              # as the map is laid out: [row from top, col]
  m <- dem_to_matrix(r)
  cands <- list(
    "t(m)"             = function(m) t(m),
    "m"                = function(m) m,
    "t(m), flip rows"  = function(m) t(m)[nrow(t(m)):1, , drop = FALSE],
    "t(m), flip cols"  = function(m) t(m)[, ncol(t(m)):1, drop = FALSE],
    "t(m), rot 180"    = function(m) t(m)[nrow(t(m)):1, ncol(t(m)):1, drop = FALSE],
    "m, flip rows"     = function(m) m[nrow(m):1, , drop = FALSE],
    "m, flip cols"     = function(m) m[, ncol(m):1, drop = FALSE],
    "m, rot 180"       = function(m) m[nrow(m):1, ncol(m):1, drop = FALSE])
  for (nm in names(cands)) {
    b <- cands[[nm]](m)
    if (identical(dim(b), dim(A)) && all(b == A)) {
      message("  matrix -> raster orientation detected: ", nm)
      assign(".MAT_TO_MAP", cands[[nm]], envir = globalenv())
      return(invisible(nm))
    }
  }
  stop("Could not determine rayshader's matrix orientation.", call. = FALSE)
}

matrix_to_rast <- function(m, dem_rast) {
  if (!exists(".MAT_TO_MAP", envir = globalenv())) detect_matrix_orientation()
  r <- terra::rast(get(".MAT_TO_MAP", envir = globalenv())(m))
  terra::ext(r) <- terra::ext(dem_rast)
  terra::crs(r) <- terra::crs(dem_rast)   # UTM 55S
  r
}

# ── rayshader OUTPUT layout, measured ─────────────────────────────────────────
# [NEW]
#   [Claude] FOUND ON THE CB TEST RUN (05_diagnostics.R, Part C): ray_shade
#   returned a 40x60 matrix for a 60x40 input, i.e. its OUTPUT is laid out
#   differently from its INPUT. Taking the output as if it were in the input
#   layout reflected every shadow and sky-view map across the reef (loggers
#   matched sky view only after "rotate 90 + mirror", and the sun calibration
#   came out with flipped handedness). master_script.R had the same problem.
#
#   Fix: on an OFF-CENTRE pillar on a NON-SQUARE grid, find the transformation
#   that puts rayshader's output back over the terrain that produced it:
#     - ambient_shade: most enclosed ground must be at the pillar's foot;
#     - ray_shade: from every sun direction the shadow must touch the pillar.
#   The chosen transformation is applied to every output immediately (in
#   make_shadow() and make_ambient()), before masking or extraction.
.output_candidates <- list(
  "as returned"       = function(o) o,
  "transpose"         = function(o) t(o),
  "flip rows"         = function(o) o[nrow(o):1, , drop = FALSE],
  "flip cols"         = function(o) o[, ncol(o):1, drop = FALSE],
  "rot 180"           = function(o) o[nrow(o):1, ncol(o):1, drop = FALSE],
  "transpose, flip rows" = function(o) t(o)[ncol(o):1, , drop = FALSE],
  "transpose, flip cols" = function(o) t(o)[, nrow(o):1, drop = FALSE],
  "transpose, rot 180"   = function(o) t(o)[ncol(o):1, nrow(o):1, drop = FALSE])

detect_output_orientation <- function(shade_fun = rayshader::ray_shade,
                                      amb_fun = rayshader::ambient_shade) {
  r <- terra::rast(nrows = 40, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 40,
                   vals = 0, crs = UTM_CRS)
  r[6:10, 8:12] <- 5                                         # pillar near the NW corner
  m   <- dem_to_matrix(r)
  gv  <- terra::values(r)[, 1] == 0
  pil <- colMeans(terra::xyFromCell(r, which(!gv)))
  xyg <- terra::xyFromCell(r, which(gv))
  dist_g <- sqrt((xyg[, 1] - pil[1])^2 + (xyg[, 2] - pil[2])^2)
  fits <- function(o) Filter(function(nm) identical(dim(.output_candidates[[nm]](o)), dim(m)),
                             names(.output_candidates))

  # ambient: most enclosed ground at the pillar's foot
  a  <- amb_fun(m, sunbreaks = 24, maxsearch = 20, zscale = 1, multicore = FALSE, progbar = FALSE)
  sa <- vapply(fits(a), function(nm) {
    v <- terra::values(matrix_to_rast(.output_candidates[[nm]](a), r))[gv, 1]
    suppressWarnings(stats::cor(dist_g, v, method = "spearman"))
  }, numeric(1))
  best_a <- names(which.max(sa))
  message(sprintf("  ambient_shade output: dims %s for input %s; using '%s' (sky view vs distance from pillar %.2f)",
                  paste(dim(a), collapse = "x"), paste(dim(m), collapse = "x"), best_a, max(sa)))
  if (max(sa) < 0.5) stop("Could not align ambient_shade output with the terrain.", call. = FALSE)

  # ray_shade: shadow must touch the pillar from every direction
  angles <- seq(0, 315, by = 45)
  outs <- lapply(angles, function(ang)
    shade_fun(m, sunaltitude = 30, sunangle = ang, zscale = 1, lambert = FALSE, progbar = FALSE))
  sd_ <- vapply(fits(outs[[1]]), function(nm) {
    mean(vapply(outs, function(o) {
      v <- terra::values(matrix_to_rast(.output_candidates[[nm]](o), r))[, 1]
      g <- which(gv & v < 0.5)
      if (length(g) < 3) return(99)
      sh <- colMeans(terra::xyFromCell(r, g)); sqrt(sum((sh - pil)^2))
    }, numeric(1)))
  }, numeric(1))
  best_s <- names(which.min(sd_))
  message(sprintf("  ray_shade output: dims %s for input %s; using '%s' (shadow centre on average %.1f m from pillar)",
                  paste(dim(outs[[1]]), collapse = "x"), paste(dim(m), collapse = "x"), best_s, min(sd_)))
  if (min(sd_) > 10) stop("Could not align ray_shade output with the terrain.", call. = FALSE)

  assign(".SHADE_FIX", .output_candidates[[best_s]], envir = globalenv())
  assign(".AMB_FIX",   .output_candidates[[best_a]], envir = globalenv())
  invisible(c(ray_shade = best_s, ambient_shade = best_a))
}

shade_fix <- function(o) {
  if (!exists(".SHADE_FIX", envir = globalenv()))
    stop("Run detect_output_orientation() first.", call. = FALSE)
  get(".SHADE_FIX", envir = globalenv())(o)
}
amb_fix <- function(o) {
  if (!exists(".AMB_FIX", envir = globalenv()))
    stop("Run detect_output_orientation() first.", call. = FALSE)
  get(".AMB_FIX", envir = globalenv())(o)
}

# One call that runs every orientation check, in the right order.
init_geometry_checks <- function(shade_fun = rayshader::ray_shade,
                                 amb_fun = rayshader::ambient_shade) {
  detect_matrix_orientation()
  check_roundtrip()
  detect_output_orientation(shade_fun, amb_fun)
  init_sun_calibration(shade_fun)
}

# Round-trip check on a real-shaped, asymmetric raster: raster -> matrix -> raster
# must give back exactly the same values in the same places.
check_roundtrip <- function() {
  r <- terra::rast(nrows = 40, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 40,
                   crs = UTM_CRS)
  terra::values(r) <- seq_len(terra::ncell(r))^1.3
  back <- matrix_to_rast(dem_to_matrix(r), r)
  if (!isTRUE(all.equal(terra::values(back)[, 1], terra::values(r)[, 1])))
    stop("Raster -> matrix -> raster round trip does not preserve the map.", call. = FALSE)
  message("  round trip raster -> matrix -> raster: OK (non-square, asymmetric test)")
  invisible(TRUE)
}


# ── Sun-direction calibration ─────────────────────────────────────────────────
# [NEW]
#   [Claude] WHY: the first real run showed ray_shade's `sunangle` is rotated
#   90° relative to the map once the DEM has gone raster -> matrix -> raster
#   (sun "east" cast a shadow north/south, not west). Rather than hard-code a
#   guess, we MEASURE the mapping on a test pillar, through exactly the same
#   path the pipeline uses, then invert it.
#
#   Bearings are compass degrees on the UTM map: 0 = north (+y), 90 = east (+x).
#   The measured mapping is  true_bearing = (s * sunangle + c) mod 360,
#   with s = +1 or -1 (handedness) and c an offset. To put the sun at a wanted
#   bearing z we therefore pass  sunangle = s * (z - c).

# Sun bearing implied by the shadow ray_shade casts for a given `sunangle`.
# The shadow falls AWAY from the sun, so the sun bearing is the opposite of the
# direction from the pillar to the shadow's centre.
observed_sun_bearing <- function(sunangle, shade_fun = rayshader::ray_shade) {
  r <- terra::rast(nrows = 201, ncols = 201, xmin = 0, xmax = 201,
                   ymin = 0, ymax = 201, vals = 0, crs = UTM_CRS)
  r[96:106, 96:106] <- 5                                   # 5 m pillar
  m  <- dem_to_matrix(r)
  s  <- shade_fix(shade_fun(m, sunaltitude = 30, sunangle = sunangle,
                            zscale = 1, lambert = FALSE, progbar = FALSE))
  sr <- matrix_to_rast(s, r)
  hv <- terra::values(r)[, 1]; sv <- terra::values(sr)[, 1]
  pillar <- terra::xyFromCell(r, which(hv > 0))
  shadow <- terra::xyFromCell(r, which(hv == 0 & sv < 0.5))  # ground cells only
  if (nrow(shadow) < 5) stop("No cast shadow found in the calibration test.", call. = FALSE)
  dx <- mean(shadow[, 1]) - mean(pillar[, 1])
  dy <- mean(shadow[, 2]) - mean(pillar[, 2])
  (atan2(-dx, -dy) * 180 / pi) %% 360                       # bearing of the sun
}

angle_diff <- function(a, b) ((a - b + 540) %% 360) - 180   # signed, in (-180, 180]

calibrate_sun_direction <- function(shade_fun = rayshader::ray_shade, tol = 10) {
  b0  <- observed_sun_bearing(0,  shade_fun)
  b90 <- observed_sun_bearing(90, shade_fun)
  d   <- angle_diff(b90, b0)
  s   <- if (abs(d - 90) < tol) 1 else if (abs(d + 90) < tol) -1 else
    stop(sprintf("Sun calibration: sunangle 0 -> %.0f°, 90 -> %.0f°; not a clean rotation.",
                 b0, b90), call. = FALSE)
  cal <- list(s = s, c = b0)
  message(sprintf("  ray_shade calibration: sunangle 0 -> sun at %.0f°, 90 -> %.0f°  (s = %+d, offset = %.0f°)",
                  b0, b90, s, b0))
  # Verify on bearings that were NOT used to fit the mapping.
  for (z in c(45, 135, 225, 315)) {
    got <- observed_sun_bearing(sunangle_for_bearing(z, cal), shade_fun)
    ok  <- abs(angle_diff(got, z)) < tol
    message(sprintf("  check: wanted sun at %3.0f°, shadow implies %3.0f° -> %s",
                    z, got, if (ok) "OK" else "FAIL"))
    if (!ok) stop("Sun calibration check failed; do not ray-trace.", call. = FALSE)
  }
  cal
}

sunangle_for_bearing <- function(bearing, cal) (cal$s * (bearing - cal$c)) %% 360

# Run once per stage; raytrace_angles() refuses to work until this has run,
# so no shadow can ever be cast with an uncalibrated direction.
init_sun_calibration <- function(shade_fun = rayshader::ray_shade) {
  if (exists(".SUN_CAL", envir = globalenv())) rm(".SUN_CAL", envir = globalenv())
  assign(".SUN_CAL", calibrate_sun_direction(shade_fun), envir = globalenv())
  invisible(get(".SUN_CAL", envir = globalenv()))
}


# ── make_ambient ──────────────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] settings read from config.R; progress bar off
#
#  Compute the ambient occlusion matrix ONCE per DEM. Same orientation/units as
#  ray_shade so it samples through the exact same extract path.
#  zscale must match what you pass to ray_shade (= dem$zscale, UTM metres).
#  multicore = TRUE would parallelise across cores but needs the parallel
#  backend set up beforehand; leaving FALSE for portability.
#
#  maxsearch is derived per-DEM from AMBIENT_MAXSEARCH_M (a
#  physical distance in METRES). pixel_size on this DEM is `zscale`, so the
#  conversion is simply maxsearch_m / zscale. This keeps the physical ray
#  distance the same on every site even when native pixel sizes differ.
make_ambient <- function(dem_matrix, zscale) {
  maxsearch_px <- max(1, round(AMBIENT_MAXSEARCH_M / zscale))
  message("  Computing ambient occlusion (sunbreaks=", AMBIENT_SUNBREAKS,
          ", maxsearch=", maxsearch_px, " px = ", AMBIENT_MAXSEARCH_M, " m) ...")
  # [Claude] amb_fix(): put rayshader's output back in the input layout
  amb_fix(rayshader::ambient_shade(heightmap = dem_matrix, sunbreaks = AMBIENT_SUNBREAKS,
                                   maxsearch = maxsearch_px, zscale = zscale,
                                   multicore = FALSE, progbar = FALSE))
}


# ── prepare_dem ───────────────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] takes a sitekey ("CB19") instead of site + year; path from config.R
#   - [Claude] skips aggregate() when the factor is 1
#   - [Claude] returns sitekey; status messages shortened
#
# Load DEM, reproject to UTM, downsample, convert to a matrix for rayshading.
# DEMs are originally in WGS84 geographic degrees (EPSG:4326).
# Reprojecting to UTM 55S puts x, y, z all in metres so zscale arithmetic is honest.
# ALSO computes the ambient occlusion matrix once and ships it in the list, so
# every hourly call reuses the same precomputed ambient layer instead of
# recomputing (ambient_shade is the slow step in this script).
#
# THREE NORMALISATION STEPS (matter for cross-site comparability):
#   1. Reef disc is a FIXED physical radius (REEF_RADIUS_M metres), not a
#      data-driven percentile. Every site gets the same physical disc area.
#   2. Off-disc pixels are pushed 100 m below reef_min BEFORE ambient_shade
#      runs, so the off-reef region can never influence ambient or direct
#      shadows on the reef. (Different photogrammetry batches left different
#      shapes of off-reef fill — rotated squares, rectangles, tight clips.
#      Without the deep-fill, ambient_shade would see those shapes as terrain.)
#   3. ray_shade / ambient_shade maxsearch are specified in METRES (constants
#      in config.R) and converted to pixels per-DEM, so ray distances are
#      physically consistent across sites.
# ── standardize_dem_crs ───────────────────────────────────────────────────────
# [NEW] Claude. Every DEM ends up in UTM_CRS (EPSG:32755) whatever it came in.
#   1. If config.R's DEM_CRS_OVERRIDE gives a CRS for this event (or its year,
#      e.g. "2017"), that CRS is ASSIGNED to the file first. Use this only when
#      the file's own CRS is missing or known to be wrong.
#   2. A DEM with no CRS and no override stops with a message.
#   3. The DEM is reprojected to UTM_CRS (as before for every site).
#   4. Checks after reprojection, each stopping with a message if it fails:
#      centre within LIZARD_BBOX (lon/lat), pixel size 0.5 mm to 5 cm, and
#      elevation range under DEM_MAX_RELIEF_M (catches cm/mm z units).
`%or%` <- function(a, b) if (is.null(a)) b else a
standardize_dem_crs <- function(dem, sitekey) {
  yr <- paste0("20", substr(sitekey, nchar(sitekey) - 1, nchar(sitekey)))
  override <- DEM_CRS_OVERRIDE[[sitekey]] %or% DEM_CRS_OVERRIDE[[yr]]
  if (!is.null(override)) {
    terra::crs(dem) <- override
    message("  CRS assigned from config.R (DEM_CRS_OVERRIDE): ", override)
  } else if (terra::crs(dem) == "") {
    stop("DEM for ", sitekey, " has no coordinate system. Find out which one it ",
         "was exported in and add it to DEM_CRS_OVERRIDE in config.R.", call. = FALSE)
  }
  src <- terra::crs(dem, describe = TRUE)
  message(sprintf("  source CRS: %s%s", src$name,
                  if (!is.na(src$code)) paste0(" (", src$authority, ":", src$code, ")") else ""))

  dem <- terra::project(dem, UTM_CRS)

  e  <- terra::ext(dem)
  ll <- sf::st_coordinates(sf::st_transform(sf::st_sfc(
          sf::st_point(c((e$xmin + e$xmax) / 2, (e$ymin + e$ymax) / 2)), crs = UTM_CRS), 4326))
  px <- ((e$xmax - e$xmin) / ncol(dem) + (e$ymax - e$ymin) / nrow(dem)) / 2
  z  <- terra::global(dem, c("min", "max"), na.rm = TRUE)
  relief <- z$max - z$min
  message(sprintf("  after reprojection: centre %.5f E, %.5f S; pixel %.4f m; relief %.2f m",
                  ll[1], -ll[2], px, relief))
  if (ll[1] < LIZARD_BBOX["lon_min"] || ll[1] > LIZARD_BBOX["lon_max"] ||
      ll[2] < LIZARD_BBOX["lat_min"] || ll[2] > LIZARD_BBOX["lat_max"])
    stop("DEM for ", sitekey, " lands at ", round(ll[1], 4), ", ", round(ll[2], 4),
         " after reprojection, not at Lizard Island. Its CRS is wrong or missing; ",
         "set the right one in DEM_CRS_OVERRIDE in config.R.", call. = FALSE)
  if (px < 0.0005 || px > 0.05)
    stop("DEM for ", sitekey, " has a pixel size of ", signif(px, 3), " m after ",
         "reprojection; expected 0.5 mm to 5 cm. Check its CRS units.", call. = FALSE)
  if (relief > DEM_MAX_RELIEF_M)
    stop("DEM for ", sitekey, " spans ", round(relief, 1), " m in elevation; expected ",
         "under ", DEM_MAX_RELIEF_M, " m. Are its heights in cm or mm?", call. = FALSE)
  dem
}

prepare_dem <- function(sitekey) {
  path <- need_file(file.path(DEM_DIR, paste0(sitekey, ".tif")))
  message("  Loading DEM: ", path)
  dem <- terra::trim(terra::rast(path))

  # reproject to UTM zone 55S so x, y, z are all in metres
  # this means zscale is meaningful (1 horizontal metre = 1 vertical metre)
  # [Claude] standardize_dem_crs(): handles DEMs delivered in other
  #   projections (e.g. the 2017 batch), and stops if the result is not at
  #   Lizard Island or not in metres.
  dem <- standardize_dem_crs(dem, sitekey)

  # ── CROP TO REEF DISC BOUNDING BOX (before aggregate, for speed) ──────────
  # The raw DEM extent can be much wider than the 6.5 m reef (e.g. CB19 covers
  # ~23 m square). After deep-fill those off-reef pixels contribute nothing,
  # but rayshader still has to march rays across them. Cropping to the disc
  # bbox + a small margin shrinks the matrix without changing physics.
  # data centroid on the native raster (mean x/y of non-NA cells)
  xy <- terra::xyFromCell(dem, which(!is.na(terra::values(dem))))
  if (nrow(xy) == 0) stop("prepare_dem: DEM has no non-NA pixels after reproject for ", sitekey)
  cx <- mean(xy[, 1]); cy <- mean(xy[, 2])
  half <- REEF_RADIUS_M + CROP_MARGIN_M
  dem  <- terra::crop(dem, terra::intersect(                 # stay inside data
    terra::ext(cx - half, cx + half, cy - half, cy + half), terra::ext(dem)))

  # ── native pixel size on the reprojected (un-aggregated) DEM ─────────────
  e <- terra::ext(dem)
  native_px <- ((e$xmax - e$xmin) / ncol(dem) + (e$ymax - e$ymin) / nrow(dem)) / 2

  # ── per-DEM downsample factor — targets REEF DISC resolution, not matrix ──
  # We want the 6.5 m reef to occupy ~reef_radius_px pixels on the downsampled
  # grid, so the disc has the same pixel density on every site. Solve for the
  # aggregation factor:
  #     target_px_size = REEF_RADIUS_M / reef_radius_px       (m per downsampled px)
  #     fact           = target_px_size / native_px_size
  # Floor at MIN_DOWNSAMPLE so a DEM already at higher resolution than the
  # target keeps its native pixels (we never upsample).
  fact <- max(MIN_DOWNSAMPLE,
              round((REEF_RADIUS_M / REEF_RADIUS_PX_TARGET) / native_px))
  dem_small <- if (fact > 1) terra::aggregate(dem, fact = fact) else dem

  # zscale = downsampled pixel size in metres (since x, y, z are all UTM metres)
  zscale <- as.numeric(native_px * fact)
  message(sprintf("  native px %.4f m, downsample x%d, zscale %.4f -> reef disc ~%d px radius",
                  native_px, fact, zscale, round(REEF_RADIUS_M / zscale)))

  # ──────────────────────────────────────────────────────────────────────────
  #  REEF MASK — fixed physical radius in METRES (replaces the 95th-percentile
  #  data-driven radius that was here before).
  #
  #  WHY FIXED, NOT DATA-DRIVEN
  #  The previous adaptive radius was sized to "where the DEM happens to have
  #  non-NA pixels", which varies by photogrammetry batch (year 17 vs 18 vs
  #  19). Wider DEMs (CB19, RS19) got bigger discs; tightly-clipped DEMs
  #  (TB18, CB18) got smaller ones. Different physical reef areas across
  #  sites = unfair cross-site comparison. Switching to a fixed REEF_RADIUS_M
  #  gives every site the SAME 13 m disc, regardless of how the DEM was clipped.
  #
  #  WHAT THIS DOES
  #  1. Find where the DEM ACTUALLY has data (data_present = !is.na).
  #  2. Centre the disc on the DATA centroid, not the bbox centre — robust to
  #     off-centre surveys (GT19 etc.).
  #  3. Convert REEF_RADIUS_M (metres) to a pixel radius for THIS DEM. zscale
  #     IS the pixel size in m on dem_small, so rad_px = REEF_RADIUS_M / zscale.
  #  4. Intersect the disc with data_present. KEY LINE: a pixel is "inside"
  #     only if it is BOTH within the disc AND a real data pixel.
  # ──────────────────────────────────────────────────────────────────────────
  nr <- nrow(dem_small); nc <- ncol(dem_small)
  # values() is ROW-MAJOR; build a row-major logical of where real data sits,
  # then reshape to an [nr x nc] matrix matching row()/col() indexing.
  v0 <- terra::values(dem_small)[, 1]
  data_present <- matrix(!is.na(v0), nrow = nr, ncol = nc, byrow = TRUE)
  rr <- row(data_present); cc <- col(data_present)
  # data centroid (mean row/col of real pixels) — robust to off-centre surveys
  row_c <- mean(rr[data_present]); col_c <- mean(cc[data_present])
  # KEY INTERSECTION: inside == within the disc AND genuinely real data.
  inside <- (sqrt((rr - row_c)^2 + (cc - col_c)^2) <= REEF_RADIUS_M / zscale) & data_present
  if (!any(inside)) stop("prepare_dem: fixed-metres disc + data_present intersection is empty for ",
                         sitekey, ". Check REEF_RADIUS_M and DEM extent.")

  # write back: everything outside the intersected mask -> NA (row-major vector)
  v <- v0; v[!as.vector(t(inside))] <- NA
  terra::values(dem_small) <- v

  # ──────────────────────────────────────────────────────────────────────────
  #  DEEP-FILL BEFORE RAYTRACING (replaces the old reef_min − relief fill)
  #  ambient_shade fires a fan of rays from each pixel and tests what fraction
  #  of the sky those rays can reach without hitting terrain. If off-reef
  #  pixels are at heights comparable to reef_min, the rays would see them as
  #  terrain and report phantom occlusion at the disc edge.
  #  Setting the off-disc fill to reef_min − 100 m guarantees those pixels are
  #  far enough below the reef that no ambient ray from inside the disc can
  #  reach them within a sensible vertical search.
  #  Crucially, this happens BEFORE make_ambient() is called below.
  # ──────────────────────────────────────────────────────────────────────────
  reef_min <- min(v, na.rm = TRUE)
  dem_small[is.na(dem_small)] <- reef_min - 100   # 100 m below; well outside any sensible search

  # matrix derived LAST, from the filled raster (rayshader needs real numbers everywhere)
  dem_matrix <- dem_to_matrix(dem_small)

  # [Claude] reef mask converted through the SAME raster -> matrix path as the
  # DEM, so it lines up with the shadow/ambient matrices whatever rayshader's
  # orientation is (it was t(inside), which assumed an orientation).
  mask_r <- dem_small
  terra::values(mask_r) <- as.numeric(as.vector(t(inside)))   # row-major, like v above
  reef_mask <- dem_to_matrix(mask_r) > 0.5

  #returns: a matrix for rayshader, a raster for coordinate extraction,
  #the zscale, the reef mask, and the ambient occlusion matrix.
  list(sitekey   = sitekey,
       matrix    = dem_matrix,
       rast      = dem_small,
       zscale    = zscale,
       # keep the reef mask so we can blank the border to NA AFTER raytracing.
       reef_mask = reef_mask,
       ambient   = make_ambient(dem_matrix, zscale))
}


# ── save_dem / load_dem ───────────────────────────────────────────────────────
# [NEW]
#   [Claude] SpatRasters can't be saved with saveRDS directly; wrap/unwrap them.
#   Replaces your in-memory dem_store, so figures can reuse the prepared DEMs
#   from a separate R session.
save_dem <- function(dem, path) {
  dem$rast <- terra::wrap(dem$rast)
  saveRDS(dem, path)
}
load_dem <- function(path) {
  dem <- readRDS(path)
  dem$rast <- terra::unwrap(dem$rast)
  dem
}


# ── load_hobo_locations ───────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] filters on sitekey; returns points already in UTM (your code
#     transformed to UTM inside each extract function)
#
# build a sf point object from geom_new coordinates for a given site/year.
# coords are in the local transverse mercator projection used when the
# geom script extracted them from the DEMs. We need to know that projection
# explicitly so st_transform can convert to UTM later.
TMERC_CRS <- paste0("+proj=tmerc +lat_0=-14.6989736557007", #this defines the coordinate system the numbers are in
                    " +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0",
                    " +datum=WGS84 +units=m +no_defs")

load_hobo_locations <- function(geom, key) {
  locs <- geom %>%
    filter(sitekey == key, !is.na(lon), !is.na(lat)) %>%
    distinct(hoboID, .keep_all = TRUE)
  if (nrow(locs) == 0) stop("No logger locations found in geom for ", key)
  message("  Loaded ", nrow(locs), " logger locations from geom for ", key)
  sf::st_transform(sf::st_as_sf(locs, coords = c("lon", "lat"), crs = TMERC_CRS),
                   UTM_CRS)
}


# ── get_site_centroid ─────────────────────────────────────────────────────────
# [YOURS]
# finds the center point of the DEM in WGS84 lat/lon degrees (EPSG:4326)
# even though the DEM is now in UTM, this reprojects back to WGS84 for suncalc
# suncalc requires geographic lat/lon — UTM coordinates are meaningless to it
get_site_centroid <- function(dem_rast) {
  e <- terra::ext(terra::project(dem_rast, "EPSG:4326"))
  list(lat = (e$ymin + e$ymax) / 2, lon = (e$xmin + e$xmax) / 2)
}


# ── sun_at_midhour ────────────────────────────────────────────────────────────
# [YOURS, EDITED] — your get_sun()
#   - [Claude] evaluated at the MIDDLE of the hour. surface_par and light_mol
#     are hourly means over hh:00–hh:59; master_script.R used the sun at hh:00.
#
#this uses a timestamp, lat and long to get the sun altitude (height above horizon in radians)
#and azimuth, which is the compass direction in radians
#suncalc requires WGS84 lat/lon — get_site_centroid() reprojects from UTM for this
sun_at_midhour <- function(hour_start, lat, lon) {
  s <- suncalc::getSunlightPosition(date = hour_start + 30 * 60, lat = lat, lon = lon)
  list(altitude = s$altitude, azimuth = s$azimuth)
}
#this determines the angle / strength of the shadows


# ── raytrace_angles ───────────────────────────────────────────────────────────
# [NEW]
#   [Claude] FIX 1 (azimuth): suncalc azimuth is 0 = SOUTH, increasing towards
#   WEST, so the compass bearing of the sun is azimuth + 180°. That bearing is
#   then converted to ray_shade's `sunangle` using the MEASURED calibration
#   (see calibrate_sun_direction above), because ray_shade's documented
#   convention did not hold on our matrices. master_script.R passed suncalc's
#   value straight in, so the thesis shadows were cast from the wrong direction.
#   [Claude] FIX 2 (refraction): the in-water sun is steeper than the in-air sun
#   (Snell's law, n = N_WATER). At 20° in air the sun is ~45° up underwater.
raytrace_angles <- function(sun, refract = TRUE) {
  if (!exists(".SUN_CAL", envir = globalenv()))
    stop("Run init_sun_calibration() before ray-tracing.", call. = FALSE)
  alt <- sun$altitude
  if (refract) alt <- pi / 2 - asin(cos(alt) / N_WATER)
  bearing <- (sun$azimuth * 180 / pi + 180) %% 360          # compass bearing of the sun
  list(alt_deg = alt * 180 / pi,
       bearing = bearing,
       az_deg  = sunangle_for_bearing(bearing, get(".SUN_CAL", envir = globalenv())))
}


# ── make_shadow ───────────────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] angles from raytrace_angles() (azimuth fix + refraction)
#   - [Claude] lambert = FALSE: loggers are weighted to sit horizontal, so we
#     want sun visibility (0–1), not a cosine weighting on the coral slope
#   - [Claude] the reef mask is applied here (was in process_timestep)
#   - [Claude] the "sun below horizon" skip moved to 02 (hours below
#     SUN_ALT_MIN_DEG are never ray-traced)
#
#this runs the raytracing in the background, sending the dem_matrix that I made with prepare_dem
#and sun positioning to do raytracing using rayshader, this returns a grid of values
#0 is fully shaded, and 1 is fully lit.
#
#  maxsearch is derived per-DEM from MAXSEARCH_M (a physical
#  distance in METRES). pixel_size on this DEM is `zscale`, so the conversion
#  is simply MAXSEARCH_M / zscale. This keeps direct-beam ray distances the
#  same on every site regardless of native pixel size.
make_shadow <- function(dem, sun, refract = TRUE) {
  a <- raytrace_angles(sun, refract)
  # [Claude] shade_fix(): put rayshader's output back in the input layout
  s <- shade_fix(rayshader::ray_shade(dem$matrix,
                            sunaltitude = a$alt_deg,
                            sunangle    = a$az_deg,
                            zscale      = dem$zscale,
                            maxsearch   = max(1, round(MAXSEARCH_M / dem$zscale)),
                            lambert     = FALSE,
                            progbar     = FALSE))
  # NA-mask the off-reef border for sampling/plotting.
  s[!dem$reef_mask] <- NA
  s
}


# ── extract_at_loggers ────────────────────────────────────────────────────────
# [YOURS, EDITED] — your extract_hobo_shadow() and extract_hobo_ambient(),
#   which were identical apart from the layer, merged into one function
#   - [Claude] buffer passed as an argument, so 2/6/10 cm can be extracted
#
#this takes the shadow grid from make_shadow, and puts it in to a spatial raster with the same coordinates as the raster, than samples each loggers exact location
#to get the logger's shadow fraction for that hour
#  Small buffer averages over a tiny disc around the logger point so we don't
#  get fooled by a single noisy pixel.
extract_at_loggers <- function(layer_matrix, dem_rast, hobo_utm, buffer_m) {
  r   <- matrix_to_rast(layer_matrix, dem_rast)
  buf <- sf::st_buffer(hobo_utm, dist = buffer_m)
  # extract mean shadow value within each buffer
  terra::extract(r, terra::vect(buf), fun = mean, na.rm = TRUE)[, 2]
}
