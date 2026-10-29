# ══════════════════════════════════════════════════════════════════════════════
#  functions_inputs.R — tides, capture-tide correction, Kd lookups, hourly table
#  Used by scripts/01_inputs.R.
#
#  PROVENANCE TAGS (every function has one):
#    [YOURS]          your code from master_script.R, logic unchanged
#    [YOURS, EDITED]  your code, with the listed changes
#    [NEW]            written by Claude
#  Comments starting "#" with no "[Claude]" marker are yours, copied from
#  master_script.R. Comments marked "[Claude]" were added by Claude.
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
#  CAPTURE-TIDE DEPTH CORRECTION
#
#  WHY
#  ───
#  dpHB / DEM heights are referenced to tide = 0 (chart datum), but the
#  photogrammetry was CAPTURED at a real, non-zero tide. The DEM therefore
#  carries the capture-day tide baked into its surface, so the computed water
#  column is too deep by exactly that tide height. The old hardcoded "- 2" only
#  ever fit ONE site (TB18); each site/year was surveyed at a different tide and
#  needs its own value.
#
#  FIX
#  ───
#     water_depth_m = (tide_cm / 100) - dpHB - capture_tide_m
#                                              ^^^^^^^^^^^^^^^^  replaces "- 2"
#
#  capture_tide_m is the tide at each site's photogrammetry timestamp, predicted
#  from the SAME tide-gauge series the predictions use, so the correction is
#  internally consistent (one tide reference throughout).
# ══════════════════════════════════════════════════════════════════════════════


# ── fit_tide_models_from_csv ──────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] directory and constituents come from config.R
#   - [Claude] ftide() called as TideHarmonics::ftide()
#
#  INTERPOLATE TIDE AT EACH CAPTURE TIME (harmonic fit from tide CSVs)
#  REPLACES an older approx()/env_long clamp. That version snapped to the nearest
#  in-window tide reading, which was unanchored noise: ALL capture dates fall
#  outside their loggers' deployment windows. Instead we fit a harmonic tide
#  model (TideHarmonics hc37) per year from the tide-gauge CSVs and predict the
#  tide at each exact capture instant.
#
#  VERIFIED on these CSVs:
#    • columns: date,time,cm,status,datetime
#    • 2017 `cm` is in METRES (0.35–2.39); 2018/2019 are CENTIMETRES (17–242,
#      9–253). Units are auto-detected and normalised to cm (replaces the old
#      hard-coded 2017 *100 patch in read_environment.R).
#    • Each CSV spans ~2–4 weeks; in-window fit is excellent (resid sd ~4 cm).
#  EXTRAPOLATION:
#    • 5 captures in-window (trustworthy); HS19/L119 2–4 d out (acceptable);
#      CB18 16.6 d out (>1 spring–neap cycle) -> WEAK, the function warns. Treat
#      CB18 as best-available or supply one external tide value for it.
fit_tide_models_from_csv <- function(dir = TIDE_CSV_DIR,
                                     years = c(2017, 2018, 2019),
                                     constituents = TIDE_CONSTITS) {
  # hc37/hc7/hc4 are LAZY-LOADED DATA objects in TideHarmonics (data/hc37.rda),
  # not exported functions — get("hc37") in the global env fails. Pull from the
  # package namespace; fall back to data() loading if needed.
  # [Claude] a vector of constituent names is passed straight to ftide()
  hc <- if (length(constituents) > 1) constituents else tryCatch(
    getExportedValue("TideHarmonics", constituents),
    error = function(e) tryCatch(
      get(constituents, envir = asNamespace("TideHarmonics")),
      error = function(e2) {
        e3 <- new.env()
        utils::data(list = constituents, package = "TideHarmonics", envir = e3)
        get(constituents, envir = e3)
      }))
  models <- list(); windows <- list()
  for (y in years) {
    f <- file.path(dir, sprintf("tides_Lizard%d.csv", y))
    if (!file.exists(f)) { warning("tide CSV not found: ", f); next }
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    d$dt <- as.POSIXct(d$datetime, format = "%d/%m/%Y %H:%M", tz = LOCAL_TZ)
    d <- d[!is.na(d$cm) & !is.na(d$dt), ]
    # unit auto-detect -> cm (Lizard range ~2.5 m = ~250 cm, so <10 means metres)
    to_cm <- if (max(abs(d$cm), na.rm = TRUE) < 10) 100 else 1
    models[[as.character(y)]]  <- TideHarmonics::ftide(d$cm * to_cm, d$dt, hc)
    windows[[as.character(y)]] <- range(d$dt)
    message(sprintf("  fitted %s tide model for %d on %d pts (%s -> %s) | units x%d -> cm",
                    paste(constituents, collapse = ","), y, nrow(d), format(min(d$dt), "%Y-%m-%d"),
                    format(max(d$dt), "%Y-%m-%d"), to_cm))
  }
  list(models = models, windows = windows)
}


# ── interp_capture_tide — harmonic prediction at the capture instant ──────────
# [YOURS, EDITED]
#   - [Claude] keyed on sitekey ("CB18") instead of site + year
#   - [Claude] dropped the unused `env` argument
#  Returns capture tide in METRES on the same datum as tide_cm.
interp_capture_tide <- function(sitekey, capture_dt, fit) {
  # 18L or 2018L -> "2018"   [Claude: taken from the last two characters of the key]
  yk <- as.character(2000L + as.integer(substr(sitekey, nchar(sitekey) - 1, nchar(sitekey))))
  m  <- fit$models[[yk]]
  if (is.null(m)) return(NA_real_)
  win <- fit$windows[[yk]]

  # extrapolation check — warn loudly if we're predicting outside the fit window
  days_out <- if (capture_dt < win[1]) as.numeric(difftime(win[1], capture_dt, units = "days"))
              else if (capture_dt > win[2]) as.numeric(difftime(capture_dt, win[2], units = "days"))
              else 0
  if (days_out > EXTRAP_WARN_DAYS)
    warning(sprintf("interp_capture_tide: %s capture is %.1f days OUTSIDE the %s ",
                    sitekey, days_out, yk),
            "tide window — harmonic prediction is extrapolated past a spring-neap ",
            "cycle and may be unreliable. Treat as best-available.")

  # predict over a ~4-min bracket and take the median (tide is flat over minutes,
  # so this is robust to predict()'s `by` units).
  pr <- tryCatch(stats::predict(m, from = capture_dt - 120, to = capture_dt + 120,
                                by = 1 / 60),
                 error = function(e) NA_real_)
  if (length(pr) == 0 || all(is.na(pr))) return(NA_real_)
  stats::median(as.numeric(pr), na.rm = TRUE) / 100   # already cm (normalised at fit) -> metres
}


# ── build_capture_tide — one capture_tide_m per site/year ─────────────────────
# [YOURS, EDITED]
#   - [Claude] your `capture_tide <- photogrammetry %>% rowwise() ...` block and
#     the CAPTURE_TIDE_OVERRIDE loop, combined into one function
#   - [Claude] records whether each value is a harmonic fit or a manual override
#   - [Claude] STOPS if a site has no capture tide; master_script.R fell back
#     to the old -2 m inside run_pipeline
#
#  MANUAL CAPTURE-TIDE OVERRIDES
#  The harmonic fit is only valid INSIDE (or within a few days of) each year's
#  CSV window. Captures far outside it extrapolate to nonsense:
#    • CB18 (Nov 02) is 16.6 d before the 2018 window -> harmonic gave ~80 m.
#    • ALL 2017 captures (Nov 07–08) are before the 2017 window (starts Nov 11)
#      -> harmonic gave 4–6 m or negative.
#  For these, paste a predicted tide (in METRES, chart datum, same datum as
#  tide_cm) from an external Lizard Island tide source at that exact timestamp.
build_capture_tide <- function(fit) {
  # [NEW] Claude, 28 Oct 2026: new 2016-aligned DEMs -> one value per SITE
  if (identical(DEPTH_REFERENCE, "survey_2016")) {
    site <- substr(SITES, 1, nchar(SITES) - 2)
    miss <- setdiff(unique(site), names(CAPTURE_TIDE_2016))
    if (length(miss) > 0)
      stop("No CAPTURE_TIDE_2016 value for site(s): ", paste(miss, collapse = ", "), call. = FALSE)
    ct <- data.frame(sitekey = SITES,
                     capture_dt = as.POSIXct(unname(CAPTURE_TIME_2016[site]),
                                             format = "%Y-%m-%d %H:%M", tz = LOCAL_TZ),
                     capture_tide_m = unname(CAPTURE_TIDE_2016[site]),
                     source = "2016 survey tide, per site (config.R CAPTURE_TIDE_2016)",
                     stringsAsFactors = FALSE)
    message("  Depth reference: 2016 survey surface, one capture tide per site")
    return(ct)
  }
  ct <- PHOTOGRAMMETRY
  ct$capture_dt <- as.POSIXct(ct$capture_dt, format = "%Y-%m-%d %H:%M", tz = LOCAL_TZ)
  ct$capture_tide_m <- mapply(interp_capture_tide, ct$sitekey, ct$capture_dt,
                              MoreArgs = list(fit = fit))
  ct$source <- "harmonic fit"
  for (k in names(CAPTURE_TIDE_OVERRIDE)) {
    hit <- ct$sitekey == k
    if (!any(hit)) stop("CAPTURE_TIDE_OVERRIDE has abbrev '", k, "' not found in PHOTOGRAMMETRY.")
    message(sprintf("  OVERRIDE %s capture tide: %.3f m -> %.3f m (manual)",
                    k, ct$capture_tide_m[hit], CAPTURE_TIDE_OVERRIDE[[k]]))
    ct$capture_tide_m[hit] <- CAPTURE_TIDE_OVERRIDE[[k]]
    ct$source[hit] <- "manual override (config.R)"
  }
  # [Claude] no silent fallback: a missing value stops the run
  missing <- setdiff(SITES, ct$sitekey[is.finite(ct$capture_tide_m)])
  if (length(missing) > 0)
    stop("No capture tide for: ", paste(missing, collapse = ", "),
         ". Add a capture time to PHOTOGRAMMETRY or a value to CAPTURE_TIDE_OVERRIDE.",
         call. = FALSE)
  ct
}


# ── read_kd_satellite — per-site KdPAR, satellite-derived ─────────────────────
# [YOURS, EDITED]
#   - [Claude] from your kdpar_lookup / KDPAR_BY_KEY block
#   - [Claude] STOPS if a site has no satellite Kd; master_script.R silently
#     used KD_CONSTANT instead
#   - [Claude] removed the KDPAR_OVERRIDE hook (hand-set per-site values)
read_kd_satellite <- function() {
  lk <- utils::read.csv(need_file(KDPAR_LOOKUP_CSV), stringsAsFactors = FALSE)
  # build a named vector keyed "SITEyy" -> chosen KdPAR, e.g. "CB18" -> 0.2358
  out <- setNames(suppressWarnings(as.numeric(lk[[KDPAR_COL]])),
                  make_key(lk$site, lk$year))
  # [Claude, 20 Oct 2026] was a stop. 2017 has no clear-sky satellite days
  #   (n_clear_days = 0), so its satellite Kd is NA. The satellite models
  #   (M2s/M3s/M4s) are then scored only on events that have one; the main
  #   models use the reference-logger Kd and are unaffected.
  missing <- setdiff(SITES, names(out)[is.finite(out)])
  if (length(missing) > 0)
    warning("No satellite ", KDPAR_COL, " for: ", paste(missing, collapse = ", "),
            ". Satellite models (M2s/M3s/M4s) will skip these events.", call. = FALSE)
  setNames(out[SITES], SITES)
}


# ── read_kd_bands — per-site 10-band Kd (spectral supplement only) ────────────
# [YOURS, EDITED]
#   - [Claude] from your KD_BANDS_BY_KEY block; returned as a long table
# Built from kd_raw_per_wavelength.csv (same file kdpar_full.R writes).
read_kd_bands <- function() {
  kb <- utils::read.csv(need_file(KD_BANDS_CSV), stringsAsFactors = FALSE)
  kb$kd <- suppressWarnings(as.numeric(kb$kd))
  kb[kb$sitekey %in% SITES, c("sitekey", "wave_nm", "kd")]
}


# ── read_env_long ─────────────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] PAR recomputed from light_lux with LUX_TO_PAR (config.R), so the
#     conversion is in this repository rather than only in read_environment.R
#   - [Claude] drops the back-calculated Kd column (your note said to remove it)
#
#env_long datasheet, and converts timestamps into a date-time object in the correct timezone for lizard island
#this sheet has kd back-calculated from the sensors, the code doesnt use it, I need to remove this
read_env_long <- function() {
  env <- utils::read.csv(need_file(ENV_LONG_CSV), stringsAsFactors = FALSE)
  env$datetime_pos <- as.POSIXct(env$datetime_pos, format = "%m/%d/%y %H:%M",
                                 tz = LOCAL_TZ)
  if (anyNA(env$datetime_pos)) stop("Some env_long timestamps failed to parse.")
  # [Claude] lux -> PAR here; falls back to the upstream light_mol if lux is absent
  if ("light_lux" %in% names(env)) {
    env$light_par <- env$light_lux * LUX_TO_PAR
  } else {
    warning("env_long has no light_lux column; using upstream light_mol as PAR.")
    env$light_par <- env$light_mol
  }
  env$Kd <- NULL   # [Claude] sensor-back-calculated Kd: never used (circular)

  # [Claude, 20 Oct 2026] a few loggers have no dpHB in env_long but do in
  #   geom_new.csv. The two files agree exactly wherever both have a value
  #   (checked on 901 loggers), so the missing ones are filled from geom.
  g  <- utils::read.csv(need_file(GEOM_CSV), stringsAsFactors = FALSE)
  g$dpHB <- suppressWarnings(as.numeric(g$dpHB))
  g  <- g[is.finite(g$dpHB), ]
  dp <- setNames(g$dpHB, paste(g$rec, g$year, g$new_unit))[!duplicated(paste(g$rec, g$year, g$new_unit))]
  env$dpHB <- suppressWarnings(as.numeric(env$dpHB))
  key  <- paste(env$site, env$year, env$unit)
  fill <- !is.finite(env$dpHB) & key %in% names(dp)
  if (any(fill)) {
    env$dpHB[fill] <- dp[key[fill]]
    message("  dpHB filled from geom_new.csv for ", length(unique(key[fill])), " logger(s): ",
            paste(unique(key[fill]), collapse = ", "))
  }

  # [NEW] Claude, 28 Oct 2026: logger heights from the NEW DEMs (config DPHB_CSV).
  #   Loggers not in the file (no position) keep the old value shifted by their
  #   event's median new-old difference, so they stay on the same reference.
  if (!is.null(DPHB_CSV) && !is.na(DPHB_CSV)) {
    nd <- utils::read.csv(need_file(DPHB_CSV), stringsAsFactors = FALSE)
    nv <- setNames(nd$dpHB, paste(nd$site, nd$year, nd$unit))
    has <- key %in% names(nv)
    old <- env$dpHB
    env$dpHB[has] <- nv[key[has]]
    ev  <- paste(env$site, env$year)
    shift <- tapply(env$dpHB[has] - old[has], ev[has], median, na.rm = TRUE)
    rest <- !has & is.finite(old) & ev %in% names(shift)
    env$dpHB[rest] <- old[rest] + shift[ev[rest]]
    message(sprintf("  dpHB from %s for %d logger(s); %d logger(s) not in it kept the old value + their event's median shift",
                    basename(DPHB_CSV), length(unique(key[has])), length(unique(key[rest]))))
  }
  env
}


# ── build_hourly_inputs ───────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] your "step2 build hourly schedule" block from run_pipeline(),
#     run once for all sites instead of inside the site loop
#   - [Claude] light_mol averaged here, over the same hour as the inputs (your
#     validation join did this separately)
#   - [Claude] sun-altitude filter moved to 02_raytrace.R (needs the DEM centroid)
#   - [Claude] reference Kd moved to 03_models.R
#
#the raw data is 5 minute intervals, but we only raytrace 1 hour to save computing.
#the section here collapses 5 min data to one row per logger per hour
build_hourly_inputs <- function(env, capture_tide, kd_sat) {
  h <- env %>%
    mutate(sitekey = make_key(site, year)) %>%
    #keep only rows belonging to this site year
    filter(sitekey %in% SITES) %>%
    #snap every 5 min timestamp to its hour
    #11:05, 11:10, ..., 11:55 all become 11:00, which makes rows share the same datetime
    mutate(datetime_pos = lubridate::floor_date(datetime_pos, TIMESTEP_INTERVAL),
           # [Claude, 20 Oct 2026] was as.character(hoboID). In 2017 the light
           #   files name loggers by site (GT01, L101) while geom_new.csv uses
           #   H01; env_long's `unit` column is the position ID in every year
           #   (identical to hoboID in 2018/2019), so it is used for the join.
           hoboID_env   = as.character(unit)) %>%
    #collapse the 5-min readings to the HOURLY MEAN per logger, so the
    #predictor (surface_par, tide, depth inputs) is mean-of-hour — matching
    #the validation join, which averages light_mol over the same hour.
    group_by(sitekey, site, year, hoboID_env, datetime_pos) %>%
    summarise(surface_par = mean(surface_par, na.rm = TRUE),
              tide_cm     = mean(tide_cm,     na.rm = TRUE),
              dpHB        = mean(dpHB,        na.rm = TRUE),
              light_mol   = mean(light_par,   na.rm = TRUE),
              n_readings  = dplyr::n(),
              .groups = "drop") %>%
    left_join(capture_tide %>% select(sitekey, capture_tide_m), by = "sitekey") %>%
    # tide-corrected water column above sensor
    # tide_cm/100 = tide above datum (m); dpHB = sensor bed elevation (m, negative down)
    # cap_tide_m  = tide at this site's photogrammetry capture time,
    #               predicted by harmonic fit (replaces the old flat -2).
    mutate(water_depth_m = tide_cm / 100 - dpHB - capture_tide_m,
           kd_sat   = unname(kd_sat[sitekey]),
           kd_const = KD_CONSTANT)
  absent <- setdiff(SITES, unique(h$sitekey))
  if (length(absent) > 0) stop("No env_long rows for: ", paste(absent, collapse = ", "))
  h
}


# ── read_geom ─────────────────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] adds a sitekey column; offsets table now lives in config.R
#
# geom_new.csv provides correct lon/lat (in local transmercator projection) and
# new_unit (= hoboID) for every logger at every site/year.
#  replaces shapefile names, which didnt have names that aligned and were never correctly joined to env_long.
read_geom <- function() {
  g <- utils::read.csv(need_file(GEOM_CSV), stringsAsFactors = FALSE)
  g <- g %>%
    distinct(rec, year, new_unit, lon, lat) %>%     # one row per logger per site/year
    rename(site = rec, hoboID = new_unit) %>%       #renames columns to match with env_long
    mutate(year = as.integer(year))
  # apply to lon/lat (geom is in local tmerc, units = metres, so dx/dy add directly)
  g %>%
    left_join(LOGGER_OFFSETS, by = c("site", "year", "hoboID")) %>%
    mutate(lon = lon + ifelse(is.na(dx_m), 0, dx_m),
           lat = lat + ifelse(is.na(dy_m), 0, dy_m)) %>%
    select(-dx_m, -dy_m) %>%
    mutate(sitekey = make_key(site, year))
}
