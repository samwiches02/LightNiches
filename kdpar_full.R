# ══════════════════════════════════════════════════════════════════════════════
#  kdpar_full.R
#  (1) Pull GBR spectral Kd(λ) (NESP TWQ 5.3, Magno-Canto et al. 2019), collapse
#      to a single photon-weighted KdPAR per site/year, write the lookup CSVs.
#  (2) Plot daily KdPAR variability across the deployment month (site-filtered).
#  (3) Map the ~1 km Kd grid around Lizard Island with site centroids labelled.
#
#  AGGREGATION (changed): we compute daily integrated KdPAR for each clear-sky
#  day in the deployment window, apply a 3-MAD outlier filter ACROSS DAYS
#  (Leys et al. 2013) with a 1.6x-median ratio backstop, then take the median.
#  This correctly drops cloud-edge contaminated days that affect all
#  wavelengths simultaneously — the prior per-wavelength pooling could not
#  catch these because the contamination isn't a per-wavelength outlier,
#  it's a per-day outlier.
#
#  Also exposed: kdpar_490corr (Morel 1988 conversion from Kd_490 to KdPAR,
#  bypasses bottom-reflectance contamination of per-band integration over
#  shallow reefs).
#
#  Validation of the underlying product against GBR in-situ observations:
#    Magno-Canto, McKinna, Robson & Fabricius (2019) Optics Express 27(20):A1350
#    Canto et al. (2021) Marine Pollution Bulletin 169:112539
#    (methods: Lee et al. 2005; McKinna et al. 2015)
#  Outlier filter:
#    Leys, Ley, Klein, Bernard & Licata (2013) J Exp Soc Psychol 49(4):764-766
#
#  REQUIREMENTS:  install.packages(c("ncdf4", "dplyr", "sf"))
# ══════════════════════════════════════════════════════════════════════════════

library(ncdf4)
library(dplyr)
library(sf)

# ── 1. CONFIG ─────────────────────────────────────────────────────────────────
USE_AGG  <- FALSE
DAP_BASE <- "https://thredds.eatlas.org.au/thredds/dodsC/NESP-TWQ-5-3_Benthic-light/kd-data"

LOCAL_DIR <- "/Users/samueleszenyi/Desktop/MEM MSc/1thesis/kd"

BOX_HALF <- 0

WAVES   <- c(412, 443, 469, 488, 531, 547, 555, 645, 667, 678)
PAR_MIN <- 400
PAR_MAX <- 700

# Kd averaging mode: "week" | "window" | "month" | "year"
# NOTE: with the new daily-first aggregation, "month" gives the MAD filter
# enough clear-sky days to engage robustly. "week" may have only 2-3 clear
# days, which is below the filter's n>=3 floor.
AVG_MODE <- "month"
PAD_DAYS <- 2
WATER_SEARCH_CELLS <- 8

# ── Outlier filter for daily-KdPAR pooling ────────────────────────────────────
# Two-stage filter applied to the vector of daily photon-weighted KdPAR
# values within the deployment window.
#   Stage 1 (MAD): 3 * MAD threshold (Leys et al. 2013). Robust but needs
#     enough spread among clean days to engage — fails when 4 of 5 days
#     are nearly identical (MAD ~= 0).
#   Stage 2 (ratio backstop): drops any single day > RATIO_THRESHOLD x
#     median. Catches the [0.21, 0.21, 0.37] case where MAD = 0.
# Set either to Inf to disable that stage.
MAD_THRESHOLD   <- 3
MAD_MIN_N       <- 3
RATIO_THRESHOLD <- 1.6

# ── Morel KdPAR conversion ────────────────────────────────────────────────────
KW_PAR     <- 0.10
CHI        <- 1.0
KD490_BAND <- 488

GEOM_PATH <- "geom_new.csv"
ENV_PATH  <- "env_long_depth_kd_longlat.csv"
SHP_DIR   <- "/Users/samueleszenyi/Desktop/MEM MSc/1thesis/shps"

# ── PLOT SITE FILTER (sections 6 & 7) ─────────────────────────────────────────
PLOT_SITES <- c("CB19", "GT19", "HS19", "N319", "TB19", "TM19")
MAP_SITES  <- "all"

# ── MAP CONFIG (section 7) ─────────────────────────────────────────────────────
MAP_WAVE <- 488
MAP_YEAR <- 2019
MAP_DATE <- as.Date("2019-11-01")
MAP_WIN  <- 10
PAD_DEG  <- 0.05

# ── 2. BUILD SITE / YEAR / DATE / COORDINATE TABLE ────────────────────────────
geom <- read.csv(GEOM_PATH, stringsAsFactors = FALSE)

site_centroids <- geom %>%
  group_by(rec, year) %>%
  summarise(local_x_m = mean(lon, na.rm = TRUE),
            local_y_m = mean(lat, na.rm = TRUE),
            n_points  = n(), .groups = "drop") %>%
  mutate(year4   = ifelse(year < 100, 2000 + year, year),
         sitekey = paste0(rec, sprintf("%02d", year %% 100)))

env <- read.csv(ENV_PATH, stringsAsFactors = FALSE)
tcol <- intersect(c("datetime_pos", "datetime", "date", "time"), names(env))[1]
if (is.na(tcol)) stop("No timestamp column found in ", ENV_PATH)

parse_dt <- function(x) {
  out <- suppressWarnings(as.POSIXct(x, format = "%m/%d/%y %H:%M", tz = "Australia/Brisbane"))
  bad <- is.na(out)
  if (any(bad)) out[bad] <- suppressWarnings(as.POSIXct(x[bad], tz = "Australia/Brisbane"))
  out
}
env$.dt <- parse_dt(env[[tcol]])

deploy_starts <- env %>%
  filter(!is.na(.dt)) %>%
  mutate(sitekey = paste0(site, sprintf("%02d", as.integer(year) %% 100))) %>%
  group_by(sitekey) %>%
  summarise(start_date = as.Date(min(.dt)),
            last_date  = as.Date(max(.dt)),
            n_obs      = n(), .groups = "drop")

read_site_coords <- function(shp_dir, sitekeys) {
  out <- lapply(sitekeys, function(key) {
    f <- file.path(shp_dir, paste0(key, ".shp"))
    if (!file.exists(f)) return(data.frame(sitekey = key, wgs_lon = NA, wgs_lat = NA))
    g  <- suppressWarnings(sf::st_read(f, quiet = TRUE))
    g4 <- sf::st_transform(g, 4326)
    ct <- sf::st_coordinates(sf::st_centroid(sf::st_union(g4)))
    data.frame(sitekey = key, wgs_lon = ct[1], wgs_lat = ct[2])
  })
  do.call(rbind, out)
}
site_coords <- read_site_coords(SHP_DIR, unique(deploy_starts$sitekey))
message("\n── Real WGS84 centroids from shapefiles ──")
print(site_coords, row.names = FALSE)

deploy_starts <- deploy_starts %>% left_join(site_coords, by = "sitekey")

sites <- site_centroids %>%
  left_join(deploy_starts, by = "sitekey") %>%
  mutate(coords_are_placeholder = is.na(wgs_lon) | is.na(wgs_lat),
         date_source = ifelse(!is.na(start_date),
                              "sensor_deployment_start", "MISSING_in_env_file")) %>%
  transmute(sitekey, rec, year = year4,
            lon = wgs_lon, lat = wgs_lat,
            local_x_m, local_y_m, n_points,
            start_date, last_date, date_source, coords_are_placeholder) %>%
  arrange(sitekey)

message("\n── Final site table (", nrow(sites), " site/years) ──")
print(as.data.frame(sites), row.names = FALSE)

# ── 3. HELPERS ────────────────────────────────────────────────────────────────
kd_url <- function(wave, year) {
  if (!is.null(LOCAL_DIR) && nzchar(LOCAL_DIR)) {
    file.path(LOCAL_DIR, sprintf("xr_kd%d_daily_%d.nc", wave, year))
  } else if (USE_AGG) {
    sprintf("%s/agg/kd%d", DAP_BASE, wave)
  } else {
    sprintf("%s/orig/kd%d/xr_kd%d_daily_%d.nc", DAP_BASE, wave, wave, year)
  }
}

find_kd_var <- function(nc, wave) {
  candidates <- c(sprintf("Kd_%d", wave), sprintf("kd_%d", wave),
                  sprintf("Kd%d", wave),  sprintf("kd%d", wave),
                  "Kd", "kd", "Kd_lee", sprintf("Kd_%d_lee", wave))
  hit <- candidates[candidates %in% names(nc$var)]
  if (length(hit) == 0) hit <- names(nc$var)[1]
  hit[1]
}

decode_dates <- function(nc, timen, year = NULL) {
  tvals <- nc$dim[[timen]]$vals
  tunit <- ncatt_get(nc, timen, "units")$value
  m <- regmatches(tunit, regexec("([A-Za-z]+)\\s+since\\s+(\\d{4}-\\d{1,2}-\\d{1,2})", tunit))[[1]]
  d_since <- as.Date(NA)
  if (length(m) == 3) {
    mult <- switch(tolower(m[2]), days=1, hours=1/24, minutes=1/1440, seconds=1/86400, 1)
    d_since <- as.Date(m[3]) + as.numeric(tvals) * mult
  }
  d_ymd <- suppressWarnings(as.Date(as.character(as.integer(tvals)), format = "%Y%m%d"))
  if (is.null(year)) return(if (all(is.na(d_since))) d_ymd else d_since)
  yr_ok <- function(d) sum(!is.na(d) & as.integer(format(d, "%Y")) == year, na.rm = TRUE)
  if (yr_ok(d_since) >= yr_ok(d_ymd)) d_since else d_ymd
}

# Two-stage outlier filter applied across DAYS (not wavelengths).
# Returns a logical KEEP vector of the same length as v, so the caller can
# apply the same mask to companion vectors (e.g. dates, other-method KdPAR).
outlier_keep <- function(v, label = "") {
  if (length(v) < 2) return(rep(TRUE, length(v)))
  keep <- rep(TRUE, length(v))
  med <- median(v)
  
  # Stage 1: MAD
  if (is.finite(MAD_THRESHOLD) && sum(keep) >= MAD_MIN_N) {
    mad_v <- mad(v[keep], constant = 1.4826)
    if (mad_v > 0) {
      mad_keep <- abs(v - med) <= MAD_THRESHOLD * mad_v
      n_drop <- sum(keep & !mad_keep)
      if (n_drop > 0) {
        message(sprintf("    [outlier MAD%s] dropped %d/%d > %g MAD (med=%.3f, mad=%.3f)",
                        ifelse(nzchar(label), paste0(" ", label), ""),
                        n_drop, length(v), MAD_THRESHOLD, med, mad_v))
        keep <- keep & mad_keep
        med <- median(v[keep])
      }
    }
  }
  
  # Stage 2: ratio backstop
  if (is.finite(RATIO_THRESHOLD) && sum(keep) >= 2 && med > 0) {
    ratio_keep <- v <= RATIO_THRESHOLD * med
    n_drop <- sum(keep & !ratio_keep)
    if (n_drop > 0) {
      message(sprintf("    [outlier RATIO%s] dropped %d/%d > %gx median %.3f",
                      ifelse(nzchar(label), paste0(" ", label), ""),
                      n_drop, length(v), RATIO_THRESHOLD, med))
      keep <- keep & ratio_keep
    }
  }
  
  keep
}

# pull nearest (or box) Kd for one wavelength, MEDIAN over a DATE RANGE.
# (Per-wavelength pooling — used internally by kdpar_one_day to fetch a
# single-day value. Daily aggregation happens in the main loop, not here.)
get_kd_one <- function(wave, year, lon, lat, date_from, date_to, pad_days,
                       diagnose = FALSE) {
  url <- kd_url(wave, year)
  nc  <- tryCatch(nc_open(url), error = function(e) NULL)
  if (is.null(nc)) { warning(sprintf("could not open %s", url))
    return(c(val = NA, ndays = 0, dcell = NA, cell_lon = NA, cell_lat = NA)) }
  on.exit(nc_close(nc), add = TRUE)
  
  dn   <- names(nc$dim)
  lonn <- dn[grep("^lon|longitude", dn, ignore.case = TRUE)][1]
  latn <- dn[grep("^lat|latitude",  dn, ignore.case = TRUE)][1]
  timen<- dn[grep("^time|date",      dn, ignore.case = TRUE)][1]
  
  lons <- nc$dim[[lonn]]$vals
  lats <- nc$dim[[latn]]$vals
  ix   <- which.min(abs(lons - lon))
  iy   <- which.min(abs(lats - lat))
  
  dates <- decode_dates(nc, timen, year)
  
  if (diagnose) {
    message(sprintf("    [diag] decoded %s..%s | pixel lon=%.3f lat=%.3f | querying %s..%s",
                    as.character(min(dates,na.rm=TRUE)), as.character(max(dates,na.rm=TRUE)),
                    lons[ix], lats[iy], as.character(date_from), as.character(date_to)))
  }
  
  lo <- date_from - pad_days; hi <- date_to + pad_days
  keep <- which(dates >= lo & dates <= hi)
  if (length(keep) == 0) return(c(val = NA, ndays = 0, dcell = NA, cell_lon = NA, cell_lat = NA))
  
  varname <- find_kd_var(nc, wave)
  dims <- nc$var[[varname]]$dim
  ord  <- sapply(dims, function(d) d$name)
  
  fillv <- tryCatch(ncatt_get(nc, varname, "_FillValue")$value, error=function(e) NA)
  sf    <- tryCatch(ncatt_get(nc, varname, "scale_factor")$value, error=function(e) 1)
  if (length(sf)==0 || is.na(sf)) sf <- 1
  fill_scaled <- if (!is.na(fillv)) fillv * sf else NA
  
  keep_finite <- function(v) {
    v <- as.numeric(v)
    ok <- is.finite(v) & v > 0 & v < 6
    if (!is.na(fill_scaled)) ok <- ok & abs(v - fill_scaled) > 1e-6
    v[ok]
  }
  
  read_box <- function(cx, cy, half) {
    x0 <- max(1, cx - half); x1 <- min(length(lons), cx + half)
    y0 <- max(1, cy - half); y1 <- min(length(lats), cy + half)
    start <- count <- integer(length(ord))
    for (k in seq_along(ord)) {
      if (ord[k] == lonn)      { start[k] <- x0;        count[k] <- x1 - x0 + 1 }
      else if (ord[k] == latn) { start[k] <- y0;        count[k] <- y1 - y0 + 1 }
      else if (ord[k] == timen){ start[k] <- min(keep); count[k] <- diff(range(keep)) + 1 }
      else                     { start[k] <- 1;         count[k] <- dims[[k]]$len }
    }
    keep_finite(ncvar_get(nc, varname, start = start, count = count))
  }
  
  box <- read_box(ix, iy, BOX_HALF)
  if (length(box) > 0)
    return(c(val = median(box), ndays = length(box), dcell = 0,
             cell_lon = lons[ix], cell_lat = lats[iy]))
  
  nlon <- length(lons); nlat <- length(lats)
  for (r in 1:WATER_SEARCH_CELLS) {
    best <- NULL; best_d <- Inf; best_lon <- NA; best_lat <- NA
    for (dx in -r:r) for (dy in -r:r) {
      if (max(abs(dx), abs(dy)) != r) next
      cx <- ix + dx; cy <- iy + dy
      if (cx < 1 || cx > nlon || cy < 1 || cy > nlat) next
      v <- read_box(cx, cy, BOX_HALF)
      if (length(v) > 0) {
        dkm <- sqrt(((lons[cx]-lons[ix]) * 111 * cos(lats[iy]*pi/180))^2 +
                      ((lats[cy]-lats[iy]) * 111)^2)
        if (dkm < best_d) { best_d <- dkm; best <- v; best_lon <- lons[cx]; best_lat <- lats[cy] }
      }
    }
    if (!is.null(best))
      return(c(val = median(best), ndays = length(best), dcell = best_d,
               cell_lon = best_lon, cell_lat = best_lat))
  }
  c(val = NA, ndays = 0, dcell = NA, cell_lon = NA, cell_lat = NA)
}

integrate_kdpar <- function(waves, kd_vals, par_min = PAR_MIN, par_max = PAR_MAX,
                            photon_weight = NULL) {
  ok <- is.finite(kd_vals)
  if (sum(ok) < 2) return(c(photon = NA, mean = NA, geom = NA, nbands = sum(ok)))
  w <- waves[ok]; k <- kd_vals[ok]
  grid <- seq(par_min, par_max, by = 1)
  k_grid <- approx(w, k, xout = grid, rule = 2)$y
  weight <- if (is.null(photon_weight)) grid else photon_weight(grid)
  weight <- weight / sum(weight)
  c(photon = sum(k_grid * weight), mean = mean(k), geom = exp(mean(log(k))), nbands = sum(ok))
}

avg_range <- function(s_row) {
  yr <- as.integer(s_row$year)
  if (AVG_MODE == "week") {
    list(from = s_row$start_date - 7, to = s_row$start_date, pad = 0)
  } else if (AVG_MODE == "month") {
    f <- as.Date(format(s_row$start_date, "%Y-%m-01"))
    list(from = f, to = seq(f, by = "month", length.out = 2)[2] - 1, pad = 0)
  } else if (AVG_MODE == "year") {
    list(from = as.Date(sprintf("%d-01-01", yr)), to = as.Date(sprintf("%d-12-31", yr)), pad = 0)
  } else {
    list(from = s_row$start_date - PAD_DAYS, to = s_row$last_date + PAD_DAYS, pad = PAD_DAYS)
  }
}

# Compute integrated KdPAR (all three weighting variants + raw Kd_490) for
# ONE specific day at ONE site. Returns NA fields if any required band is
# missing for that day.
kdpar_one_day <- function(s_row, the_date) {
  yr <- as.integer(s_row$year)
  kd_vals <- sapply(WAVES, function(w)
    unname(get_kd_one(w, yr, s_row$lon, s_row$lat, the_date, the_date, 0,
                      diagnose = FALSE)["val"]))
  integ <- integrate_kdpar(WAVES, kd_vals)
  list(
    photon = unname(integ["photon"]),
    mean   = unname(integ["mean"]),
    geom   = unname(integ["geom"]),
    kd_490 = if (KD490_BAND %in% WAVES) kd_vals[which(WAVES == KD490_BAND)] else NA_real_
  )
}

# ── 4. MAIN LOOP — one KdPAR per site/year ────────────────────────────────────
# NEW FLOW: aggregate by DAY first, then filter outliers across days, then
# take median across surviving days. This is the order that catches
# cloud-edge contaminated days (the HS19 spike) because the contamination
# inflates ALL wavelengths on that day, making it a per-day outlier even
# when it isn't a per-wavelength outlier.
results <- list(); raw_long <- list(); daily_long <- list()

for (i in seq_len(nrow(sites))) {
  s <- sites[i, ]
  if (is.na(s$start_date)) { message(sprintf("[%s] no date — skipping", s$sitekey)); next }
  if (is.na(s$lon) || is.na(s$lat)) { message(sprintf("[%s] no coords — skipping", s$sitekey)); next }
  yr <- as.integer(s$year)
  rng <- avg_range(s)
  
  message(sprintf("[%s] avg=%s  %s .. %s  (%.4f, %.4f)  year=%d",
                  s$sitekey, AVG_MODE, rng$from, rng$to, s$lon, s$lat, yr))
  
  # ── 4a. DAILY AGGREGATION ───────────────────────────────────────────────────
  window_days <- seq(rng$from, rng$to, by = "day")
  message(sprintf("   [daily] querying %d days in window...", length(window_days)))
  
  daily_photon <- numeric(0); daily_mean <- numeric(0)
  daily_geom   <- numeric(0); daily_kd490 <- numeric(0)
  daily_dates  <- as.Date(character(0))
  for (d in window_days) {
    d <- as.Date(d, origin = "1970-01-01")
    one <- kdpar_one_day(s, d)
    if (is.finite(one$photon) && is.finite(one$mean) && is.finite(one$geom)) {
      daily_photon <- c(daily_photon, one$photon)
      daily_mean   <- c(daily_mean,   one$mean)
      daily_geom   <- c(daily_geom,   one$geom)
      daily_kd490  <- c(daily_kd490,  one$kd_490)
      daily_dates  <- c(daily_dates,  d)
    }
  }
  message(sprintf("   [daily] %d/%d days returned valid KdPAR",
                  length(daily_photon), length(window_days)))
  
  # ── 4b. OUTLIER FILTER ACROSS DAYS ──────────────────────────────────────────
  # Run the filter on the photon-weighted series (the most physically
  # meaningful aggregation) and apply the same day-mask to all four series
  # so they remain internally consistent.
  keep <- outlier_keep(daily_photon, label = sprintf("[daily %s]", s$sitekey))
  n_before <- length(daily_photon)
  daily_photon <- daily_photon[keep]
  daily_mean   <- daily_mean[keep]
  daily_geom   <- daily_geom[keep]
  daily_kd490  <- daily_kd490[keep]
  daily_dates  <- daily_dates[keep]
  if (length(daily_photon) < n_before)
    message(sprintf("   [daily] %d/%d days survived outlier filter",
                    length(daily_photon), n_before))
  
  # ── 4c. POOLED VALUES (median across surviving days) ────────────────────────
  pooled_photon <- if (length(daily_photon) > 0) median(daily_photon) else NA_real_
  pooled_mean   <- if (length(daily_mean)   > 0) median(daily_mean)   else NA_real_
  pooled_geom   <- if (length(daily_geom)   > 0) median(daily_geom)   else NA_real_
  pooled_kd490  <- if (length(daily_kd490)  > 0) median(daily_kd490, na.rm = TRUE) else NA_real_
  pooled_490corr <- if (is.finite(pooled_kd490)) KW_PAR + CHI * pooled_kd490 else NA_real_
  
  message(sprintf("   -> KdPAR (n=%d days): photon=%.3f mean=%.3f geom=%.3f | Kd_490=%.3f -> 490corr=%.3f",
                  length(daily_photon), pooled_photon, pooled_mean, pooled_geom,
                  pooled_kd490, pooled_490corr))
  
  results[[i]] <- data.frame(
    sitekey       = s$sitekey,
    start_date    = s$start_date,
    last_date     = s$last_date,
    lon           = s$lon,
    lat           = s$lat,
    kdpar_photon  = pooled_photon,
    kdpar_mean    = pooled_mean,
    kdpar_geom    = pooled_geom,
    kd_490        = pooled_kd490,
    kdpar_490corr = pooled_490corr,
    nbands        = length(WAVES),
    n_clear_days  = length(daily_photon),
    avg_mode      = AVG_MODE,
    stringsAsFactors = FALSE
  )
  
  # per-band raw values (across the full unfiltered window) for traceability
  for (j in seq_along(WAVES)) {
    out_w <- get_kd_one(WAVES[j], yr, s$lon, s$lat, rng$from, rng$to, rng$pad)
    raw_long[[length(raw_long) + 1]] <- data.frame(
      sitekey = s$sitekey, wave_nm = WAVES[j],
      kd = unname(out_w["val"]), n_days = unname(out_w["ndays"]),
      water_cell_km = round(unname(out_w["dcell"]), 2))
  }
  
  # per-day record (post-filter) for the variability plot
  if (length(daily_photon) > 0) {
    daily_long[[length(daily_long) + 1]] <- data.frame(
      sitekey    = s$sitekey,
      date       = daily_dates,
      kdpar      = daily_photon,
      start_date = s$start_date,
      stringsAsFactors = FALSE)
  }
}

kdpar_table <- bind_rows(results)
kd_raw      <- bind_rows(raw_long)
daily       <- bind_rows(daily_long)

cat("\n================ KdPAR per site/year (", AVG_MODE,
    ", daily-first + ", MAD_THRESHOLD, "MAD + ", RATIO_THRESHOLD,
    "x ratio filter) ================\n", sep="")
print(kdpar_table, row.names = FALSE, digits = 4)

write.csv(kdpar_table, "kdpar_by_site_year.csv", row.names = FALSE)
write.csv(kd_raw,      "kd_raw_per_wavelength.csv", row.names = FALSE)
message("\nWrote kdpar_by_site_year.csv and kd_raw_per_wavelength.csv")

# ── 5. LOOKUP for master_script.R ─────────────────────────────────────────────
kdpar_lookup <- kdpar_table %>%
  mutate(site = substr(sitekey, 1, nchar(sitekey) - 2),
         year = 2000 + as.integer(substr(sitekey, nchar(sitekey) - 1,
                                         nchar(sitekey)))) %>%
  select(sitekey, site, year,
         kdpar_photon, kdpar_mean, kdpar_geom, kdpar_490corr,
         kd_490, nbands, n_clear_days)
write.csv(kdpar_lookup, "kdpar_lookup.csv", row.names = FALSE)
message("Wrote kdpar_lookup.csv (daily-first aggregation)")

cat("\n── Kd comparison across methods (m^-1) ──\n")
cmp <- kdpar_lookup %>%
  select(sitekey, n_clear_days, kdpar_photon, kdpar_mean, kdpar_geom, kdpar_490corr) %>%
  mutate(ratio_photon_to_490corr = round(kdpar_photon / kdpar_490corr, 2))
print(as.data.frame(cmp), row.names = FALSE, digits = 3)
cat("\nn_clear_days < 5 = limited sample; per-site KdPAR uncertainty is\n",
    "higher at those sites and should be noted in methods.\n", sep = "")

# ══════════════════════════════════════════════════════════════════════════════
#  6. DAILY KdPAR VARIABILITY PLOT  (post-filter — what was actually used)
# ══════════════════════════════════════════════════════════════════════════════
plot_keys_all <- sites %>% filter(!is.na(start_date), !is.na(lon)) %>% pull(sitekey)

plot_keys <- if (identical(PLOT_SITES, "all")) plot_keys_all else intersect(PLOT_SITES, plot_keys_all)
bad <- setdiff(if (identical(PLOT_SITES,"all")) character(0) else PLOT_SITES, plot_keys_all)
if (length(bad)) warning("PLOT_SITES not in run: ", paste(bad, collapse=", "))
if (length(plot_keys) == 0) stop("PLOT_SITES matched nothing.")
message("\nDaily-variability plot includes: ", paste(plot_keys, collapse = ", "))

daily_plot <- daily[daily$sitekey %in% plot_keys & is.finite(daily$kdpar), ]

if (nrow(daily_plot) > 0) {
  keys <- sort(unique(daily_plot$sitekey))
  pal  <- setNames(hcl.colors(max(length(keys),2), "Dark 3")[seq_along(keys)], keys)
  
  dev.new(width = 11, height = 6)
  par(mar = c(5, 5, 4, 8))
  daily_plot$dom <- as.integer(format(daily_plot$date, "%d"))
  xr <- range(daily_plot$dom); yrng <- range(daily_plot$kdpar, na.rm = TRUE)
  
  plot(NA, xlim = xr, ylim = yrng,
       xlab = "Day of deployment month", ylab = "Daily KdPAR (m^-1)",
       main = paste0("Daily KdPAR variability (post-filter) — ",
                     paste(keys, collapse=", "),
                     "\n(points = surviving clear-sky days; solid = median; X = deployment day)"),
       cex.main = 0.9)
  for (key in keys) {
    sub <- daily_plot[daily_plot$sitekey == key, ]; sub <- sub[order(sub$dom), ]
    points(sub$dom, sub$kdpar, col = pal[key], pch = 19, cex = 0.9)
    lines(sub$dom, sub$kdpar, col = pal[key], lwd = 1, lty = 3)
    abline(h = median(sub$kdpar, na.rm = TRUE), col = pal[key], lwd = 1.5)
    points(as.integer(format(sub$start_date[1], "%d")),
           median(sub$kdpar, na.rm = TRUE),
           col = pal[key], pch = 4, cex = 1.8, lwd = 2)
  }
  par(xpd = TRUE)
  legend(xr[2] + diff(xr) * 0.03, yrng[2], legend = keys, col = pal[keys],
         pch = 19, lwd = 1.5, bty = "n", title = "Site/year")
  par(xpd = FALSE)
}

cat("\n── KdPAR stability within the deployment window (post-filter) ──\n")
daily %>% group_by(sitekey) %>%
  summarise(n_clear_days = n(),
            kdpar_median = round(median(kdpar),3),
            kdpar_mean   = round(mean(kdpar),3),
            kdpar_sd     = round(sd(kdpar),3),
            kdpar_cv_pct = round(100*sd(kdpar)/mean(kdpar),1),
            kdpar_min    = round(min(kdpar),3),
            kdpar_max    = round(max(kdpar),3),
            .groups="drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  7. MAP — unchanged
# ══════════════════════════════════════════════════════════════════════════════
map_src <- sites[!is.na(sites$lon) & !is.na(sites$lat), c("sitekey","lon","lat")]
map_src$site <- substr(map_src$sitekey, 1, nchar(map_src$sitekey) - 2)
if (!identical(MAP_SITES, "all")) {
  want_sites <- unique(substr(MAP_SITES, 1, nchar(MAP_SITES) - 2))
  map_src <- map_src[map_src$site %in% want_sites, ]
}
site_pts <- aggregate(cbind(lon, lat) ~ site, data = map_src, FUN = mean)
cat("\nSites on map:\n"); print(site_pts, row.names = FALSE)

nc_path <- file.path(LOCAL_DIR, sprintf("xr_kd%d_daily_%d.nc", MAP_WAVE, MAP_YEAR))
nc <- nc_open(nc_path)
dn   <- names(nc$dim)
lonn <- dn[grep("^lon|longitude", dn, ignore.case = TRUE)][1]
latn <- dn[grep("^lat|latitude",  dn, ignore.case = TRUE)][1]
timen<- dn[grep("^time|date",      dn, ignore.case = TRUE)][1]
lons <- nc$dim[[lonn]]$vals; lats <- nc$dim[[latn]]$vals

xlim <- range(site_pts$lon) + c(-PAD_DEG, PAD_DEG)
ylim <- range(site_pts$lat) + c(-PAD_DEG, PAD_DEG)
xi <- which(lons >= xlim[1] & lons <= xlim[2])
yi <- which(lats >= ylim[1] & lats <= ylim[2])
lon_sub <- lons[xi]; lat_sub <- lats[yi]

dates <- decode_dates(nc, timen, MAP_YEAR)
win <- which(dates >= MAP_DATE - MAP_WIN & dates <= MAP_DATE + MAP_WIN)
if (length(win) == 0) { nc_close(nc); stop("No dates near MAP_DATE in the file.") }

vname <- find_kd_var(nc, MAP_WAVE)
dims <- nc$var[[vname]]$dim; ord <- sapply(dims, function(d) d$name)
start <- count <- integer(length(ord))
for (k in seq_along(ord)) {
  if (ord[k]==lonn)      { start[k] <- min(xi);  count[k] <- length(xi) }
  else if (ord[k]==latn) { start[k] <- min(yi);  count[k] <- length(yi) }
  else if (ord[k]==timen){ start[k] <- min(win); count[k] <- diff(range(win)) + 1 }
  else                   { start[k] <- 1;        count[k] <- dims[[k]]$len }
}
arr <- ncvar_get(nc, vname, start = start, count = count)
nc_close(nc)

arr[!is.finite(arr) | arr <= 0 | arr > 6] <- NA
if (length(dim(arr)) == 3) {
  kd <- apply(arr, c(1, 2), function(v) if (all(is.na(v))) NA else median(v, na.rm = TRUE))
} else {
  kd <- arr
}

if (lat_sub[1] > lat_sub[length(lat_sub)]) { lat_sub <- rev(lat_sub); kd <- kd[, ncol(kd):1] }
if (lon_sub[1] > lon_sub[length(lon_sub)]) { lon_sub <- rev(lon_sub); kd <- kd[nrow(kd):1, ] }

has_kd <- any(is.finite(kd))
if (!has_kd)
  warning("Map Kd field is all-NA over this window.")

nearest_idx <- function(v, x) which.min(abs(v - x))
sampled_cell <- function(px, py) {
  ix <- nearest_idx(lon_sub, px); iy <- nearest_idx(lat_sub, py)
  if (!is.na(kd[ix, iy])) return(c(lon = lon_sub[ix], lat = lat_sub[iy], dist_km = 0))
  nlon <- length(lon_sub); nlat <- length(lat_sub)
  for (r in 1:WATER_SEARCH_CELLS) {
    best <- NULL; best_d <- Inf
    for (dx in -r:r) for (dy in -r:r) {
      if (max(abs(dx), abs(dy)) != r) next
      cx <- ix + dx; cy <- iy + dy
      if (cx < 1 || cx > nlon || cy < 1 || cy > nlat) next
      if (!is.na(kd[cx, cy])) {
        dkm <- sqrt(((lon_sub[cx]-lon_sub[ix]) * 111 * cos(lat_sub[iy]*pi/180))^2 +
                      ((lat_sub[cy]-lat_sub[iy]) * 111)^2)
        if (dkm < best_d) { best_d <- dkm; best <- c(lon = lon_sub[cx], lat = lat_sub[cy], dist_km = dkm) }
      }
    }
    if (!is.null(best)) return(best)
  }
  c(lon = NA, lat = NA, dist_km = NA)
}

dev.new(width = 9, height = 8)
par(mar = c(5, 5, 4, 6))
pal <- hcl.colors(100, "YlGnBu", rev = TRUE)
zlim <- if (has_kd) range(kd, na.rm = TRUE) else c(0, 1)

image(lon_sub, lat_sub, kd, col = pal, zlim = zlim,
      xlab = "Longitude", ylab = "Latitude",
      main = paste0("GBR-Kd ~1 km grid around Lizard Island\n",
                    MAP_WAVE, " nm, ", format(MAP_DATE), " +/-", MAP_WIN, " d"),
      useRaster = FALSE, asp = 1/cos(mean(lat_sub)*pi/180))
image(lon_sub, lat_sub, ifelse(is.na(kd), 1, NA), col = "grey80", add = TRUE, useRaster = FALSE)

dx <- median(diff(lon_sub)); dy <- median(diff(lat_sub))
abline(v = c(lon_sub - dx/2, max(lon_sub) + dx/2), col = "grey90", lwd = 0.4)
abline(h = c(lat_sub - dy/2, max(lat_sub) + dy/2), col = "grey90", lwd = 0.4)

cell_log <- data.frame()
for (i in seq_len(nrow(site_pts))) {
  px <- site_pts$lon[i]; py <- site_pts$lat[i]
  sc <- sampled_cell(px, py)
  nx <- lon_sub[nearest_idx(lon_sub, px)]; ny <- lat_sub[nearest_idx(lat_sub, py)]
  if (!is.na(sc["lon"])) {
    if (abs(sc["lon"] - nx) > dx/4 || abs(sc["lat"] - ny) > dy/4) {
      rect(nx - dx/2, ny - dy/2, nx + dx/2, ny + dy/2, border = "grey50", lwd = 1, lty = 2)
      segments(nx, ny, sc["lon"], sc["lat"], col = "grey50", lwd = 1, lty = 3)
    }
    rect(sc["lon"] - dx/2, sc["lat"] - dy/2, sc["lon"] + dx/2, sc["lat"] + dy/2,
         border = "red", lwd = 2)
  }
  points(px, py, pch = 21, bg = "red", col = "white", cex = 1.4, lwd = 1.5)
  text(px, py, labels = site_pts$site[i], pos = 3, offset = 0.5, cex = 0.9, font = 2)
  cell_log <- rbind(cell_log, data.frame(site = site_pts$site[i],
                                         cell_lon = round(unname(sc["lon"]), 4),
                                         cell_lat = round(unname(sc["lat"]), 4),
                                         spiral_km = round(unname(sc["dist_km"]), 2)))
}

if (has_kd) {
  uxr <- par("usr")[1:2]; uyr <- par("usr")[3:4]
  strip_x <- uxr[2] + diff(uxr)*0.02
  sy <- seq(uyr[1] + diff(uyr)*0.1, uyr[1] + diff(uyr)*0.5, length.out = 100)
  par(xpd = TRUE)
  for (k in 1:99) rect(strip_x, sy[k], strip_x + diff(uxr)*0.015, sy[k+1], col = pal[k], border = NA)
  text(strip_x, max(sy) + diff(uyr)*0.03, paste0("Kd", MAP_WAVE), cex = 0.8, font = 2, adj = 0)
  text(strip_x + diff(uxr)*0.02, min(sy), round(zlim[1],2), cex = 0.7, adj = 0)
  text(strip_x + diff(uxr)*0.02, max(sy), round(zlim[2],2), cex = 0.7, adj = 0)
  par(xpd = FALSE)
}

cell_log$cell_id <- paste(cell_log$cell_lon, cell_log$cell_lat, sep = "_")
n_unique <- length(unique(cell_log$cell_id[!is.na(cell_log$cell_lon)]))
cat("\n── Sampled cell per site (red outline) ──\n")
print(cell_log[, c("site","cell_lon","cell_lat","spiral_km")], row.names = FALSE)
cat(sprintf("\n%d sites resolve to %d UNIQUE grid cells.\n", nrow(cell_log), n_unique))
dupe <- cell_log[duplicated(cell_log$cell_id) | duplicated(cell_log$cell_id, fromLast = TRUE), ]
dupe <- dupe[!is.na(dupe$cell_lon), ]
if (nrow(dupe) > 0) {
  cat("Sites SHARING a cell:\n")
  for (id in unique(dupe$cell_id))
    cat("  ", paste(dupe$site[dupe$cell_id == id], collapse = " + "), "\n")
} else {
  cat("All mapped sites fall in distinct cells.\n")
}

cat(sprintf("\nGrid resolution: %.3f km (lon) x %.3f km (lat)\n",
            dx*111*cos(mean(lat_sub)*pi/180), dy*111))

