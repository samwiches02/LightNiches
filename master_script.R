# ══════════════════════════════════════════════════════════════════════════════
#  LIGHT NICHES
#  Predicts underwater PAR at HOBO logger locations
# ══════════════════════════════════════════════════════════════════════════════
#kdpar_full.R ──► kdpar_lookup.csv ──┐   (per-site KdPAR, satellite-derived,
#(satellite Kd, run separately)      │    independent of the loggers)
#                                    │
#DEM .tif  ──┐                       │
#geom_new ──┼──► prepare_dem()       │   env_long ──► capture-tide correction
#env_long  ──┘   • reproject UTM 55S │   • interp tide at photogrammetry capture
#• downsample→target │   • water_depth = tide − dpHB − capture_tide
#• circle-mask reef  │
#• ambient_shade()   │
#(computed ONCE)     │
#│                            │
#▼                            ▼
#get_site_centroid()           per logger, per hour:
#  (UTM→WGS84 for suncalc)        Kd (from lookup)  +  water_depth_m
#│                            │
#┌────┴───────── for each hour ────┴──────────┐
#▼                                            ▼
#get_sun() (alt+az)                    calc_water_attenuation()
#│                                   exp(−Kd × depth)  → 0–1
#▼                                   [OR spectral per-band sum, toggle]
#make_shadow() = ray_shade() → 0–1 matrix       │
#│            ambient matrix (reused)          │
#▼                                             │
#extract_hobo_shadow() / extract_hobo_ambient()  │
#(sample both matrices at logger points)         │
#│                                             │
#▼                                             │
#eff_light = W_DIRECT·shadow + W_AMBIENT·ambient │
#[OR shadow only, if USE_AMBIENT=F]   │
#W_AMBIENT may be static OR altitude-dependent
#└──────────────────┬──────────────────────────┘
#▼
#pred_par = surface_par × eff_light × water_atten
#│
#▼
#Validation:  pred_par ↔ light_mol
#(per-logger panels, shadow maps, pooled DLI, co-fit)



# ── LIBRARIES ─────────────────────────────────────────────────────────────────

library(terra)       # raster/vector GIS package
library(raster)      # needed by rayshader
library(rayshader)   # ray_shade() casts shadow maps; ambient_shade() casts AO maps
library(rgl)         # used by rayshader
library(suncalc)     # used for get sunlight position, uses long/lat/time
library(sf)          # read shapefiles, reproject etc.
library(dplyr)       # data wrangling
library(tidyr)       # data reshaping
library(lubridate)   # date-time helper
library(readr)       # tidyverse uses this
library(TideHarmonics)  # harmonic tide fit for the capture-tide depth correction



# ── WORKING DIRECTORY ─────────────────────────────────────────────────────────
#set working directory, need to change to github repository

setwd("~/Desktop/MEM MSc/1thesis/LightNiches-main/output")


# ── LOAD INPUT DATA ───────────────────────────────────────────────────────────

#env_long datasheet, and converts timestamps into a date-time object in the correct timezone for lizard island
#this sheet has kd back-calculated from the sensors, the code doesnt use it, I need to remove this
env_long <- read.csv("env_long_depth_kd_longlat.csv") %>%
  mutate(datetime_pos = as.POSIXct(datetime_pos, format = "%m/%d/%y %H:%M", tz = "Australia/Brisbane"))


# geom_new.csv provides correct lon/lat (in local transmercator projection) and
# new_unit (= hoboID) for every logger at every site/year.
# this replaces shapefile names, which didnt have names that aligned and were never correctly joined to env_long.
geom <- read.csv("geom_new.csv") %>%
  distinct(rec, year, new_unit, lon, lat) %>%     # one row per logger per site/year
  rename(site = rec, hoboID = new_unit) %>% #renames columns to match with env_long
  mutate(year = as.integer(year))


# ══════════════════════════════════════════════════════════════════════════════
#  CAPTURE-TIDE DEPTH CORRECTION  (built here, applied inside run_pipeline)
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
#
#  NOTE ON KEYS: env_long / geom use site = "TB"/"CB"/"RS" and year = 18/19, so
#  the photogrammetry table below is keyed the SAME way. If it were keyed on the
#  long names ("turtle_beach", 2018) the left_join would silently yield NA and
#  every water_depth_m would become NA.
# ══════════════════════════════════════════════════════════════════════════════
#####
# ── PHOTOGRAMMETRY CAPTURE TIMES (local Brisbane time) ────────────────────────
photogrammetry <- tibble::tribble(
  ~site, ~year, ~capture_dt,         ~abbrev,
  "TB",  18L,   "2018-11-30 14:39",  "TB18",
  "CB",  18L,   "2018-11-02 08:32",  "CB18",
  "RS",  18L,   "2018-12-02 16:15",  "RS18",
  "TB",  19L,   "2019-11-03 13:54",  "TB19",
  "CB",  19L,   "2019-11-03 09:09",  "CB19",
  "GT",  19L,   "2019-11-08 10:07",  "GT19",
  "GT",  17L,   "2017-11-08 11:22",  "GT17",
  "HS",  17L,   "2017-11-08 16:43",  "HS17",
  "HS",  19L,   "2019-11-02 10:24",  "HS19",
  "L1",  17L,   "2017-11-08 16:08",  "L117",
  "L1",  19L,   "2019-11-03 16:20",  "L119",
  "N3",  17L,   "2017-11-07 10:01",  "N317",
  "N3",  18L,   "2018-11-28 13:09",  "N318",
  "N3",  19L,   "2019-11-03 10:24",  "N319",
  "SE",  17L,   "2017-11-07 09:13",  "SE17",
  "SE",  18L,   "2018-11-25 08:55",  "SE18",
  "TM",  18L,   "2018-11-23 09:24",  "TM18",
  "TM",  19L,   "2019-11-12 10:43",  "TM19"
) %>%
  mutate(capture_dt = as.POSIXct(capture_dt,
                                 format = "%Y-%m-%d %H:%M",
                                 tz     = "Australia/Brisbane"))

# ── INTERPOLATE TIDE AT EACH CAPTURE TIME (harmonic fit from tide CSVs) ────────
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

# Directory holding tides_Lizard{2017,2018,2019}.csv. Absolute path is fine even
# though setwd() points at /output — just make sure the 3 CSVs sit in here.
TIDE_CSV_DIR     <- "/Users/samueleszenyi/Desktop/MEM MSc/1thesis/LightNiches-main/tides"
TIDE_CONSTITS    <- "hc37"   # read_environment.R used hc37
EXTRAP_WARN_DAYS <- 7        # warn if a capture is >this far outside the fit window

# ── fit one harmonic model per year, normalising units to CENTIMETRES ─────────
fit_tide_models_from_csv <- function(dir = TIDE_CSV_DIR,
                                     years = c(2017, 2018, 2019),
                                     constituents = TIDE_CONSTITS) {
  # hc37/hc7/hc4 are LAZY-LOADED DATA objects in TideHarmonics (data/hc37.rda),
  # not exported functions — get("hc37") in the global env fails. Pull from the
  # package namespace; fall back to data() loading if needed.
  hc <- tryCatch(
    getExportedValue("TideHarmonics", constituents),
    error = function(e) tryCatch(
      get(constituents, envir = asNamespace("TideHarmonics")),
      error = function(e2) {
        e3 <- new.env()
        utils::data(list = constituents, package = "TideHarmonics", envir = e3)
        get(constituents, envir = e3)
      }
    )
  )
  models <- list(); windows <- list()
  for (y in years) {
    f <- file.path(dir, sprintf("tides_Lizard%d.csv", y))
    if (!file.exists(f)) { warning("tide CSV not found: ", f); next }
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    d$dt <- as.POSIXct(d$datetime, format = "%d/%m/%Y %H:%M", tz = "Australia/Brisbane")
    d <- d[!is.na(d$cm) & !is.na(d$dt), ]
    # unit auto-detect -> cm (Lizard range ~2.5 m = ~250 cm, so <10 means metres)
    to_cm <- if (max(abs(d$cm), na.rm = TRUE) < 10) 100 else 1
    d$cm_norm <- d$cm * to_cm
    models[[as.character(y)]]  <- ftide(d$cm_norm, d$dt, hc)
    windows[[as.character(y)]] <- range(d$dt)
    message(sprintf("  fitted %s tide model for %d on %d pts (%s -> %s) | units x%d -> cm",
                    constituents, y, nrow(d),
                    format(min(d$dt), "%Y-%m-%d"), format(max(d$dt), "%Y-%m-%d"), to_cm))
  }
  list(models = models, windows = windows)
}

message("\n── Fitting harmonic tide models from CSVs ──")
tide_fit <- fit_tide_models_from_csv()

# 18L or 2018L -> "2018"
.tide_year_key <- function(year_i) {
  yi <- as.integer(year_i); if (yi < 100) yi <- 2000L + yi; as.character(yi)
}

# ── interp_capture_tide — harmonic prediction at the capture instant ──────────
#  Returns capture tide in METRES on the same datum as tide_cm. `env` is accepted
#  but UNUSED so the existing call site (which passes env_long) needs no change.
interp_capture_tide <- function(site_i, year_i, capture_dt_i, env = NULL,
                                fit = tide_fit) {
  yk <- .tide_year_key(year_i)
  if (is.null(fit$models[[yk]])) {
    warning("interp_capture_tide: no model for year ", yk, " (site ", site_i,
            "). Returning NA — run_pipeline falls back to -2 m.")
    return(NA_real_)
  }
  m   <- fit$models[[yk]]
  win <- fit$windows[[yk]]
  
  # extrapolation check — warn loudly if we're predicting outside the fit window
  days_out <- if (capture_dt_i < win[1])
    as.numeric(difftime(win[1], capture_dt_i, units = "days"))
  else if (capture_dt_i > win[2])
    as.numeric(difftime(capture_dt_i, win[2], units = "days"))
  else 0
  if (days_out > EXTRAP_WARN_DAYS)
    warning(sprintf("interp_capture_tide: %s capture is %.1f days OUTSIDE the %s ",
                    site_i, days_out, yk),
            "tide window — harmonic prediction is extrapolated past a spring-neap ",
            "cycle and may be unreliable. Treat as best-available.")
  
  # predict over a ~4-min bracket and take the median (tide is flat over minutes,
  # so this is robust to predict()'s `by` units).
  from <- capture_dt_i - as.difftime(2, units = "mins")
  to   <- capture_dt_i + as.difftime(2, units = "mins")
  pr <- tryCatch(stats::predict(m, from = from, to = to, by = 1/60),
                 error = function(e) { warning("predict() failed for ", site_i,
                                               " ", yk, ": ", conditionMessage(e)); NA_real_ })
  if (length(pr) == 0 || all(is.na(pr))) return(NA_real_)
  
  tide_cm <- stats::median(as.numeric(pr), na.rm = TRUE)   # already cm (normalised at fit)
  tide_cm / 100                                            # -> metres
}

capture_tide <- photogrammetry %>%
  rowwise() %>%
  mutate(capture_tide_m = interp_capture_tide(site, year, capture_dt, env_long)) %>%
  ungroup()

# ── MANUAL CAPTURE-TIDE OVERRIDES ─────────────────────────────────────────────
#  The harmonic fit is only valid INSIDE (or within a few days of) each year's
#  CSV window. Captures far outside it extrapolate to nonsense:
#    • CB18 (Nov 02) is 16.6 d before the 2018 window -> harmonic gave ~80 m.
#    • ALL 2017 captures (Nov 07–08) are before the 2017 window (starts Nov 11)
#      -> harmonic gave 4–6 m or negative. (You don't run 2017 now, but these
#      would be wrong if you did.)
#  For these, paste a predicted tide (in METRES, chart datum, same datum as
#  tide_cm) from an external Lizard Island tide source at that exact timestamp.
#  Keyed by `abbrev`. Any abbrev NOT listed keeps its harmonic value.
#  Leave a value as NA to force the old -2 m fallback in run_pipeline instead.
CAPTURE_TIDE_OVERRIDE <- c(
  # "CB18" = 0.00   # <-- REQUIRED for CB18: external Lizard tide @ 2018-11-02 08:32
  # "GT17" = 0.00,  # only needed if you add 2017 sites
  # "N317" = 0.00,
  # "SE17" = 0.00,
  # "HS17" = 0.00,
  # "L117" = 0.00
)

if (length(CAPTURE_TIDE_OVERRIDE) > 0) {
  for (ab in names(CAPTURE_TIDE_OVERRIDE)) {
    hit <- capture_tide$abbrev == ab
    if (any(hit)) {
      old <- capture_tide$capture_tide_m[hit]
      capture_tide$capture_tide_m[hit] <- CAPTURE_TIDE_OVERRIDE[[ab]]
      message(sprintf("  OVERRIDE %s capture tide: %.3f m -> %.3f m (manual)",
                      ab, old, CAPTURE_TIDE_OVERRIDE[[ab]]))
    } else {
      warning("CAPTURE_TIDE_OVERRIDE has abbrev '", ab,
              "' not found in capture_tide.")
    }
  }
}

# ── LOOKUP TABLE — one capture_tide_m per site/year, subtracted in the depth
#    formula instead of the flat "- 2". Printed so you can eyeball it before runs.
message("\n── Capture-tide lookup table ──")
print(capture_tide %>% select(abbrev, site, year, capture_dt, capture_tide_m))

# ── PART C SANITY CHECK — capture tides should fall in each year's tide range ──
#  Prints each year's observed tide range (from env_long) next to the predicted
#  capture tides. In-window captures should sit inside [lo_m, hi_m]; CB18 (the
#  extrapolated one, already warned above) is the value to eyeball.
local({
  yr_range <- env_long %>% filter(!is.na(tide_cm)) %>%
    mutate(yr = as.integer(format(datetime_pos, "%Y"))) %>%
    group_by(yr) %>%
    summarise(lo_m = min(tide_cm) / 100, hi_m = max(tide_cm) / 100, .groups = "drop")
  message("\n── Part C: tide range per year (m) vs capture tides ──")
  print(as.data.frame(yr_range), row.names = FALSE)
  print(capture_tide %>% select(abbrev, capture_dt, capture_tide_m) %>% as.data.frame(),
        row.names = FALSE)
})

#####

# ── FILE-PATH ROOTS ───────────────────────────────────────────────────────────
#create filepath for DEMS
DEM_DIR <- "/Users/samueleszenyi/Desktop/MEM MSc/1thesis/dems"

# ── constants ─────────────────────────────────────────────────────────────────

KD_CONSTANT  <- 0.32   # fallback per-site, AND value used when USE_CONSTANT_KD=TRUE
USE_CONSTANT_KD <- FALSE   # TRUE = ignore the lookup; apply KD_CONSTANT to every site

# Per-site implied Kd from residuals
# Beer-Lambert: pred = surface * eff * exp(-Kd_used * z)
#               meas = surface * eff * exp(-Kd_true * z)
# => Kd_true = Kd_used + ln(pred/meas) / z

kd_implied <- validation %>%
  filter(is.finite(pred_par), is.finite(light_mol),
         pred_par > 50, light_mol > 50,         # daylight only
         water_depth_m > 1) %>%                 # avoid divide-by-tiny
  mutate(kd_inferred = Kd + log(pred_par / light_mol) / water_depth_m) %>%
  group_by(sitekey) %>%
  summarise(kd_used    = first(Kd),
            kd_median  = median(kd_inferred, na.rm = TRUE),
            kd_mean    = mean(kd_inferred,   na.rm = TRUE),
            kd_p25     = quantile(kd_inferred, 0.25, na.rm = TRUE),
            kd_p75     = quantile(kd_inferred, 0.75, na.rm = TRUE),
            n          = n(),
            .groups = "drop")

print(kd_implied)

# which KdPAR variant to apply:
#   "kdpar_photon" = photon-weighted full-band integral
#   "kdpar_mean"   = unweighted mean of the 10 bands
#   "kdpar_geom"   = geometric mean, red-outlier-resistant
KDPAR_COL <- "kdpar_photon"

kdpar_lookup <- read.csv("kdpar_lookup.csv", stringsAsFactors = FALSE)
# build a named vector keyed "SITEyy" -> chosen KdPAR, e.g. "CB18" -> 0.2358
kdpar_lookup$key <- paste0(kdpar_lookup$site, sprintf("%02d", kdpar_lookup$year %% 100))
KDPAR_BY_KEY <- setNames(kdpar_lookup[[KDPAR_COL]], kdpar_lookup$key)
if (USE_CONSTANT_KD) {
  KDPAR_BY_KEY[] <- KD_CONSTANT
  message("USE_CONSTANT_KD ON — overriding every site to Kd = ", KD_CONSTANT)
}
message("Loaded per-site KdPAR (", KDPAR_COL, ") for ", length(KDPAR_BY_KEY),
        " site/years. Fallback Kd = ", KD_CONSTANT)

# ── per-site KdPAR printout (sorted, with all three variants for comparison) ──
local({
  show <- kdpar_lookup[, c("key", "kdpar_photon", "kdpar_mean", "kdpar_geom", "nbands")]
  show <- show[order(show$key), ]
  show[["<- in use"]] <- ifelse(names(KDPAR_BY_KEY)[match(show$key, names(KDPAR_BY_KEY))] == show$key,
                                paste0("[", KDPAR_COL, "]"), "")
  cat("\n── Per-site KdPAR (active column = ", KDPAR_COL, ") ──\n", sep = "")
  print(show, row.names = FALSE, digits = 3)
})


# ── PER-SITE Kd OVERRIDE HOOK ─────────────────────────────────────────────────
#  Manually replace the satellite KdPAR for specific site/years. Useful when the
#  depth-gradient regression on a site disagrees with the satellite value (e.g.
#  CB19 wants ~0.55 but satellite gives 0.24). Keys are "SITEyy", same format
#  as KDPAR_BY_KEY. Any key NOT listed keeps its satellite value.
#
#  Skipped entirely when USE_CONSTANT_KD = TRUE (constant Kd deliberately
#  ignores everything) and when USE_REFERENCE_KD = TRUE (reference-derived
#  Kd(t) supersedes the per-site scalar). Only takes effect on satellite-Kd
#  runs and on the spectral path's fallback.
KDPAR_OVERRIDE <- c(
  "CB19" = 0.28, #mean
  "GT19" = 0.30, #mean
  "HS19" = 0.27, #median
  "N319" = 0.32, #median
  "TB19" = 0.58, #median
  "TM19" = 0.42 #median
)

if (length(KDPAR_OVERRIDE) > 0 && !USE_CONSTANT_KD && !USE_REFERENCE_KD) {
  for (k in names(KDPAR_OVERRIDE)) {
    if (k %in% names(KDPAR_BY_KEY)) {
      message(sprintf("  Kd OVERRIDE %s: %.3f -> %.3f",
                      k, KDPAR_BY_KEY[[k]], KDPAR_OVERRIDE[[k]]))
      KDPAR_BY_KEY[[k]] <- KDPAR_OVERRIDE[[k]]
    } else {
      warning("KDPAR_OVERRIDE has key '", k, "' not in KDPAR_BY_KEY — ignored.")
    }
  }
} else if (length(KDPAR_OVERRIDE) > 0) {
  message("KDPAR_OVERRIDE set but USE_CONSTANT_KD or USE_REFERENCE_KD is on — ignored.")
}

# ── SPECTRAL (depth-dependent) Kd toggle ──────────────────────────────────────
# FALSE = single photon-weighted KdPAR (your validated main result).
# TRUE  = attenuate each of the 10 bands at its own Kd and sum surviving PAR;
#         depth-dependence then emerges automatically (Canto method). Run BOTH
#         and compare the DLI plot / RMSE — this is a sensitivity analysis, not
#         a replacement, until the comparison says otherwise.
USE_SPECTRAL_KD <- FALSE

WAVES_PAR <- c(412, 443, 469, 488, 531, 547, 555, 645, 667, 678)  # must match kdpar_full.R

# Per-site 10-band Kd vector, keyed "SITEyy" -> numeric(10) in WAVES_PAR order.
# Built from kd_raw_per_wavelength.csv (same file kdpar_full.R writes).
kd_raw_tab <- read.csv("kd_raw_per_wavelength.csv", stringsAsFactors = FALSE)
KD_BANDS_BY_KEY <- local({
  out <- list()
  for (k in unique(kd_raw_tab$sitekey)) {
    sub <- kd_raw_tab[kd_raw_tab$sitekey == k, ]
    sub <- sub[match(WAVES_PAR, sub$wave_nm), ]   # enforce WAVES_PAR order
    out[[k]] <- sub$kd
  }
  out
})
if (USE_SPECTRAL_KD) {
  ok_keys <- names(KD_BANDS_BY_KEY)[sapply(KD_BANDS_BY_KEY, function(v) !any(is.na(v)))]
  message("Spectral Kd ON. Complete 10-band vectors for: ",
          paste(ok_keys, collapse = ", "))
} else {
  message("Spectral Kd OFF — using single ", KDPAR_COL, ".")
}



# ── REFERENCE-LOGGER TIME-VARYING Kd (sensitivity / alternate attenuation) ────
# FALSE = use the satellite per-site KdPAR (your main, independent result).
# TRUE  = derive Kd(t) per site by inverting Beer-Lambert on one or more open,
#         known-depth reference loggers:
#           Kd(t) = -ln(ref_light / surface_par) / ref_depth
#         The reference logger is an INPUT (like surface_par) and is EXCLUDED
#         from validation, so this is not circular — but it requires a genuinely
#         UNSHADED reference at correct depth, and ref_light & surface_par in the
#         SAME units. CHECK the printed transmittance diagnostics before trusting.
USE_REFERENCE_KD <- FALSE


# ── enforce mutual exclusion: constant Kd overrides the others ────────────────
if (USE_CONSTANT_KD) {
  USE_SPECTRAL_KD  <- FALSE
  USE_REFERENCE_KD <- FALSE
  message("USE_CONSTANT_KD ON — disabling USE_SPECTRAL_KD and USE_REFERENCE_KD")
}

# Reference loggers — 1..N candidates per site. Per-hour Kd from each candidate
# is computed independently, then aggregated across candidates by REF_AGG_FUN.
# Median is more robust than mean if one ref has unmodelled shadow at some hour
# (a single inflated Kd value gets ignored; mean would drag the aggregate up).
N_REF_LOGGERS <- 5
REF_AGG_FUN   <- function(x) stats::median(x, na.rm = TRUE)

# Paste the find_reference_loggers(validation) output here.
# REFERENCE_LOGGER must be a LIST (not a named vector) — the validation-drop
# step uses lengths()/unlist() to expand multiple IDs per site.
#
# BOOTSTRAP WORKFLOW (because the new finder requires light_mol from validation,
# and validation depends on the Kd method, the recommendation depends on the
# current run's Kd method):
#   1. Set USE_REFERENCE_KD <- FALSE for ONE run (pred_par becomes satellite-Kd).
#   2. Run the script; finder prints REFERENCE_LOGGER <- list(...) recommendation.
#   3. Paste that block in here.
#   4. Set USE_REFERENCE_KD <- TRUE and re-run.
#   5. Don't re-run the finder unless something structural changes.
REFERENCE_LOGGER <- list(
  "CB19" = c("H11", "H57"),
  "GT19" = c("H27", "H37"),
  "HS19" = c("H57", "H52"),
  "N319" = c("H44", "H01"),
  "TB19" = c("H47", "H18", "H52", "H05", "H46"),
  "TM19" = c("H67", "H28", "H47")
)


TIMESTEP_INTERVAL <- "1 hour"    # the timestep interval to measure

# ── DOWNSAMPLING: target the REEF DISC's pixel resolution, not the matrix ─────
# OLD behaviour: TARGET_PIXELS = 1000 targeted the FULL matrix width. But the
# raw DEMs cover different physical extents (CB18 ~13 m wide, CB19 ~23 m wide,
# rectangular sites even more), so 1000 px on the full matrix gave VERY
# different pixel resolutions inside the reef itself:
#   CB18: 6.5 m reef in 1058 px matrix -> reef ~ 525 px radius
#   CB19: 6.5 m reef in  961 px matrix -> reef ~ 265 px radius (HALF as fine!)
# Different reef pixel densities = different raytrace fidelity per site, which
# is the very inconsistency we are trying to eliminate.
#
# NEW behaviour: target the REEF DISC's pixel radius directly. Set
# REEF_RADIUS_PX_TARGET (px on the downsampled grid that the 6.5 m reef should
# fill), and derive the downsample factor per-DEM from the NATIVE pixel size so
# the disc ends up at that resolution on every site.
#
#   target_px_size = REEF_RADIUS_M / REEF_RADIUS_PX_TARGET   (metres per px)
#   fact           = round( target_px_size / native_px_size )
#
# At REEF_RADIUS_PX_TARGET = 500 the disc is 13 m / 1000 px = 1.3 cm per px on
# every site — fine enough to resolve individual coral heads, consistent across
# sites. Bigger value = finer disc = slower raytrace (cost scales ~N^2).
# Smaller value = coarser disc but faster. floor-guarded at MIN_DOWNSAMPLE so
# we never upsample (a low-res DEM stays at native resolution).
REEF_RADIUS_PX_TARGET <- 600   # target pixel radius of the reef disc
MIN_DOWNSAMPLE        <- 1     # never upsample; 1 = keep native if already finer
# zscale still adapts automatically — it is native_pixel_size * (the per-DEM
# factor), so coarser or finer DEMs each get the correct vertical scaling.
#
# logic (unchanged, for reference):
#a zscale of 1 will always be wrong right now, because the DEMs are always 13m wide
#BUT the native pixel resolution will change depending on the nature of the thing used to take the dem image
#so if it is 6346 pixels, 1 native pixel = 0.002m, if downsampled by 10x, each pixel is 0.02m(2cm)
#so the zscale is 0.02
#in rayshader documentation - zscale = xy pixel spaceing / z unit ratio
#both are now in meters afer UTM, the pixel size becomes meters
#for zscale to = 1, each pixel would need to represent 1 meter, so that means a 13m wide dem would have 13 pixels across
#bad
#if the resolution is as high as 2mm no matter which projection, then zscale will be a small fraction of 1

#Buffer to extract shadow and ambient values (always matching) in meters
BUFFER_M <- 0.7

# ── REEF DISC & RAY-SEARCH — ALL IN METRES (replaces pixel-based equivalents) ─
#
#  WHY THESE EXIST AS METRE-CONSTANTS NOW
#  ───────────────────────────────────────
#  The DEMs were produced in three separate photogrammetry batches (years 17,
#  18, 19), each with different clipping conventions. Some are tightly clipped
#  to the ~13 m reef (TB18, CB18); some carry a wider rotated-square fill area
#  around the reef (CB19, RS19, several other 2019 files); GT19 is a rectangle.
#  Without normalisation, every site has a different *physical* footprint
#  inside the pipeline:
#
#  • OLD REEF MASK was sized to the 95th percentile of where the DEM happened
#    to have non-NA pixels. On a wide unclipped DEM that percentile sits
#    outside the real reef, so the disc inflates to 16–18 m. On a tightly
#    clipped DEM it stays at ~13 m. Different sites = different reef areas.
#
#  • OLD MAXSEARCH was in PIXELS. ray_shade(maxsearch = 700) reaches a
#    different *physical* distance per DEM (700 px × different pixel size).
#    Combined with the variable disc, the shadow physics differed site-to-site
#    for reasons unrelated to actual reef geometry.
#
#  • OLD OFF-REEF FILL was reef_min − relief (one full relief-depth below).
#    Mostly fine but on tilted rotated-square fills the corners could still
#    sit close enough to the reef minimum to cast phantom occlusion on the
#    ambient layer near the disc edge.
#
#  NEW (METRES-BASED) BEHAVIOUR
#  ────────────────────────────
#  • REEF_RADIUS_M     : every site gets a 6.5 m radius (~13 m diameter) disc
#                        centred on the data centroid. Identical physical reef
#                        area on every site, regardless of how the DEM was
#                        processed.
#  • MAXSEARCH_M       : direct-beam rays travel the same physical distance
#                        per site. Converted to pixels per-DEM via zscale
#                        inside make_shadow().
#  • AMBIENT_MAXSEARCH_M: same for the ambient layer; converted in make_ambient().
#  • DEEP-FILL         : off-disc pixels are pushed 100 m below reef_min BEFORE
#                        ambient_shade runs, so the off-reef region cannot
#                        contribute to ambient sky-view or to direct shadows.
#                        The raytrace output is still NA-blanked outside the
#                        disc afterwards so plots show a clean circular reef.
#
#  Result: TB18 / CB18 (already tight) should look the same as before; CB19
#  and other unclipped 2019 DEMs collapse to the same clean disc.
REEF_RADIUS_M       <- 6.5   # ~13 m diameter disc; matches survey footprint
MAXSEARCH_M         <- 13  # direct-beam ray distance (m). At TB18's zscale
# (~0.0123 m/px on a 1062px DEM) this is ~700 px,
# matching the old MAXSEARCH_FACT=700 — so TB18
# reproduces previous shadows; other sites become
# consistent with it.



# ── AMBIENT OCCLUSION constants ───────────────────────────────────────────────
# ambient_shade() models DIFFUSE skylight: for each pixel it fires a fan of rays
# in a circle (sunbreaks) at a range of altitudes (anglebreaks) and measures what
# fraction of the sky it can see. A pixel in a crevice sees little sky -> low
# value (darker); an open flat or ridge sees the whole sky -> ~1. This is the
# physically-motivated "diffuse floor" so that even a fully ray_shade-shaded
# crevice is not pitch black underwater.
# It depends ONLY on geometry + zscale (NOT sun position), so it is computed ONCE
# per DEM and reused for every hour.
USE_AMBIENT <- TRUE # FALSE = pure direct-beam (eff_light = shadow only), sensitivity test

AMBIENT_MAXSEARCH_M <- 13   # ambient ray distance (m). ~500 px on TB18.
AMBIENT_SUNBREAKS  <- 30    # rays sent in a circle around each pixel (rayshader default is 24)

# maxsearch is set per-DEM inside make_ambient from AMBIENT_MAXSEARCH_M above.
# static-mode weights (used when USE_ALT_DEPENDENT_WEIGHTS = FALSE)
W_DIRECT  <- 0.3
W_AMBIENT <- 0.7   # NB: W_DIRECT + W_AMBIENT must equal 1


# per-site best static W_AMBIENT (much faster than 3-D, isolates per-site bias) (doesnt work before validation is created)
if (exists("validation")){
validation %>%
  group_by(sitekey) %>%
  group_modify(~ {
    g <- seq(0, 0.95, by = 0.05)
    rmse <- sapply(g, function(w) {
      eff <- (1-w)*.x$shadow_fraction + w*.x$ambient_fraction
      pred <- .x$surface_par * eff * exp(-.x$Kd * .x$water_depth_m)
      sqrt(mean((pred - .x$light_mol)^2, na.rm = TRUE))
    })
    data.frame(best_w = g[which.min(rmse)], rmse = min(rmse))
  }) %>% print()
}



# ── DIRECT/AMBIENT BLEND ──────────────────────────────────────────────────────
# eff_light = w_direct·shadow_fraction + w_ambient·ambient_fraction
# Two modes, controlled by USE_ALT_DEPENDENT_WEIGHTS:
#
#  FALSE: static weights (W_DIRECT, W_AMBIENT below). Simple, single co-fit
#         possible via fit_ambient_weights(). Sum must equal 1.
#
#  TRUE:  W_AMBIENT varies with sun altitude (physically motivated). The diffuse
#         fraction of incoming PAR genuinely changes with sun height — clear-sky
#         empirical and radiative-transfer evidence agree.
#
#  Functional form:
#     W_AMBIENT(alt) = W_AMBIENT_MIN + (MAX-MIN) * (1-sin(alt))^P
#         At zenith (alt=90°) -> W_AMBIENT_MIN  (clear-sky PAR diffuse fraction
#         ~0.15-0.20, Spitters et al. 1986 AgForMet 38:217).
#         At low sun -> W_AMBIENT_MAX (diffuse dominates: longer atm. path
#         per Kasten & Young 1989 ApplOpt 28:4735; stronger Rayleigh scatter
#         in PAR per Weiss & Norman 1985 AgForMet 34:205; clear-sky diffuse
#         fraction rises to 0.6-0.8 below 15° per Reindl et al. 1990 SolEnergy
#         45:1, Erbs et al. 1982 SolEnergy 28:293). Underwater the low-sun
#         direct beam has a longer refracted path (Snell, n=1.33; Kirk 2011
#         Light&Photosynthesis Ch.6; Mobley 1994 Light&Water) so it attenuates
#         faster than near-isotropic diffuse, reinforcing the shift.
#
#  WHY THE TWO MODES COEXIST: static is the simple sensitivity comparator (and
#  is what fit_ambient_weights() can search). Altitude-dependent is the
#  physically defensible default. Flip USE_ALT_DEPENDENT_WEIGHTS to compare.
USE_ALT_DEPENDENT_WEIGHTS <- TRUE   # FALSE = use static W_DIRECT/W_AMBIENT

USE_DEPTH_DEPENDENT_WEIGHTS <- TRUE   # NEW toggle, sits alongside USE_ALT_DEPENDENT_WEIGHTS



# altitude-dependent-mode endpoints (used when USE_ALT_DEPENDENT_WEIGHTS = TRUE)
W_AMBIENT_MIN <- 0.10  # diffuse floor: shallow + high sun
W_AMBIENT_MAX <- 0.80  # diffuse cap:  deep or low sun

W_AMBIENT_P   <- 2     # altitude curve shape
#   1 = linear in (1-sin α) — ramps up early
#   2 = stays near MIN until ~30°, then ramps up
# P=2 preserves spatial discrimination through midday
# (when sun is high and direct beam dominates the signal
# we actually want to separate loggers on), and only
# lets diffuse take over at low sun when the direct
# beam is genuinely a minor share of incoming PAR.

TAU_SCALE     <- 0.8  # optical-depth e-folding (τ at which depth-term = 0.63)
ALT_EXP   <- 0.8    # damp altitude — only matters at very low sun
DEPTH_EXP <- 1.2  # amplify depth — even modest depth pulls ambient up



# helper: convert sun altitude (radians) -> W_AMBIENT
# pmax/sin guards against negative altitudes (sun below horizon); pmin/pmax
# clamps to [0,1] so the result is always a valid weight.
# Combined altitude + optical-depth diffuse weight.
# Set USE_ALT_DEPENDENT_WEIGHTS=FALSE to zero out the altitude term.
# Set USE_DEPTH_DEPENDENT_WEIGHTS=FALSE to zero out the depth term.
# Both off => static W_AMBIENT (handled in process_timestep, not here).

w_ambient_combined <- function(sun_alt_rad, water_depth_m, Kd) {
  # If both axes are off, behave like static mode (return W_AMBIENT).
  # Otherwise fall through to the alt × depth blend.
  if (!USE_ALT_DEPENDENT_WEIGHTS && !USE_DEPTH_DEPENDENT_WEIGHTS) {
    return(rep(W_AMBIENT, length(water_depth_m)))
  }
  
  A <- if (USE_ALT_DEPENDENT_WEIGHTS) {
    s <- sin(pmax(sun_alt_rad, 0))
    (1 - s)^W_AMBIENT_P
  } else 0
  
  D <- if (USE_DEPTH_DEPENDENT_WEIGHTS) {
    tau <- pmax(Kd * water_depth_m, 0)
    1 - exp(-tau / TAU_SCALE)
  } else 0
  
  combined <- 1 - (1 - A)^ALT_EXP * (1 - D)^DEPTH_EXP
  w <- W_AMBIENT_MIN + (W_AMBIENT_MAX - W_AMBIENT_MIN) * combined
  pmin(pmax(w, 0), 1)
}

# print the curve once at startup so you can verify the chosen P/MIN/MAX
# produces sensible weights across the day
if (USE_AMBIENT && (USE_ALT_DEPENDENT_WEIGHTS || USE_DEPTH_DEPENDENT_WEIGHTS)) {
  message("\n── W_AMBIENT preview (alt × depth) ──")
  alts <- c(90, 60, 30, 15) * pi/180
  depths <- c(1, 2, 3, 4)
  Kd_demo <- KD_CONSTANT
  tab <- outer(alts, depths, function(a, z)
    round(w_ambient_combined(a, z, Kd_demo), 2))
  dimnames(tab) <- list(paste0("alt", round(alts*180/pi), "°"),
                        paste0("z=", depths, "m"))
  print(tab)
} else {
  message("Static weights mode: W_DIRECT=", W_DIRECT, " W_AMBIENT=", W_AMBIENT)
}

# ── plot-title weight label ───────────────────────────────────────────────────
# Returns a short string describing the current weighting scheme. Used in every
# plot title so figures self-document which mode produced them (you can't mix
# up a static-w_amb=0.6 figure with an alt-dependent figure when comparing).
weight_label <- function() {
  if (!USE_AMBIENT) return("ambient OFF (shadow only)")
  parts <- c()
  if (USE_ALT_DEPENDENT_WEIGHTS)   parts <- c(parts, "alt-dep")
  if (USE_DEPTH_DEPENDENT_WEIGHTS) parts <- c(parts, sprintf("depth-dep(τ=%.1f)", TAU_SCALE))
  if (length(parts) == 0)
    sprintf("static [w_dir=%.2f w_amb=%.2f]", W_DIRECT, W_AMBIENT)
  else
    sprintf("w=%s [%.2f→%.2f, p=%g]",
            paste(parts, collapse="+"), W_AMBIENT_MIN, W_AMBIENT_MAX, W_AMBIENT_P)
}

# ══════════════════════════════════════════════════════════════════════════════
#  TUNE ALTITUDE-DEPENDENT DIFFUSE WEIGHTING
#  3-D grid over (W_AMBIENT_MIN, W_AMBIENT_MAX, W_AMBIENT_P). Per-row Kd held at
#  its satellite value (so this isolates the curve shape from attenuation).
#  Honours USE_AMBIENT and the spectral-Kd path. Pooled across all site/years.
# ══════════════════════════════════════════════════════════════════════════════
fit_alt_weights <- function(validation,
                            min_grid = seq(0.3, 0.70, by = 0.05),
                            max_grid = seq(0.70, 0.90, by = 0.05),
                            p_grid   = c(1, 1.5, 2, 3, 4)) {
  
  if (!USE_AMBIENT) {
    warning("USE_AMBIENT is FALSE — eff_light is shadow only, weights have no ",
            "effect. Set USE_AMBIENT <- TRUE before tuning the curve.")
    return(invisible(NULL))
  }
  
  fitdat <- validation %>%
    filter(!is.na(shadow_fraction), !is.na(ambient_fraction),
           !is.na(light_mol), !is.na(surface_par),
           !is.na(water_depth_m), !is.na(sun_altitude), !is.na(Kd))
  
  message("Tuning on ", nrow(fitdat), " rows across ",
          length(unique(fitdat$sitekey)), " site/years")
  
  # Pre-compute the depth-attenuation term ONCE — it doesn't depend on weights,
  # so recomputing it inside the triple loop would be ~thousands of wasted exp().
  water_atten <- exp(-fitdat$Kd * fitdat$water_depth_m)
  sin_alt     <- sin(pmax(fitdat$sun_altitude, 0))
  
  grid <- expand.grid(w_min = min_grid, w_max = max_grid, p = p_grid)
  grid <- grid[grid$w_max > grid$w_min, ]            # enforce MAX > MIN
  grid$rmse   <- NA_real_
  grid$r2_1to1 <- NA_real_
  grid$bias   <- NA_real_
  
  for (i in seq_len(nrow(grid))) {
    w_amb <- grid$w_min[i] +
      (grid$w_max[i] - grid$w_min[i]) * (1 - sin_alt)^grid$p[i]
    w_amb <- pmin(pmax(w_amb, 0), 1)
    eff   <- (1 - w_amb) * fitdat$shadow_fraction + w_amb * fitdat$ambient_fraction
    pred  <- fitdat$surface_par * eff * water_atten
    
    res <- pred - fitdat$light_mol
    grid$rmse[i]    <- sqrt(mean(res^2, na.rm = TRUE))
    grid$bias[i]    <- mean(res, na.rm = TRUE)
    ss_res <- sum(res^2, na.rm = TRUE)
    ss_tot <- sum((fitdat$light_mol - mean(fitdat$light_mol, na.rm = TRUE))^2,
                  na.rm = TRUE)
    grid$r2_1to1[i] <- 1 - ss_res / ss_tot
  }
  
  grid <- grid[order(grid$rmse), ]
  best <- grid[1, ]
  
  cat("\n=====================================================\n")
  cat("ALTITUDE-DEPENDENT DIFFUSE-WEIGHT CO-FIT\n")
  cat("=====================================================\n")
  cat("Searched ", nrow(grid), " (MIN, MAX, P) combinations\n", sep = "")
  cat("Best W_AMBIENT_MIN :", best$w_min, "\n")
  cat("Best W_AMBIENT_MAX :", best$w_max, "\n")
  cat("Best W_AMBIENT_P   :", best$p, "\n")
  cat("Best RMSE          :", round(best$rmse, 2), "\n")
  cat("R^2 vs 1:1         :", round(best$r2_1to1, 3), "\n")
  cat("Bias               :", round(best$bias, 2),
      "(+ve = over-predicts)\n\n")
  cat("Top 10 combinations:\n")
  print(head(grid, 10), row.names = FALSE, digits = 3)
  
  # ── visualise the response surface ──────────────────────────────────────────
  # one MIN x MAX heatmap per P value, so you see which P regime dominates
  dev.new(width = 4 * length(p_grid), height = 4)
  par(mfrow = c(1, length(p_grid)), mar = c(4, 4, 3, 1))
  for (pv in p_grid) {
    sub <- grid[grid$p == pv, ]
    if (nrow(sub) < 4) { plot.new(); title(paste0("P=", pv, " (sparse)")); next }
    mat <- matrix(NA_real_, length(min_grid), length(max_grid),
                  dimnames = list(min_grid, max_grid))
    for (k in seq_len(nrow(sub)))
      mat[as.character(sub$w_min[k]), as.character(sub$w_max[k])] <- sub$rmse[k]
    image(min_grid, max_grid, mat, col = hcl.colors(100, "viridis"),
          xlab = "W_AMBIENT_MIN", ylab = "W_AMBIENT_MAX",
          main = sprintf("P=%g | RMSE", pv))
    contour(min_grid, max_grid, mat, add = TRUE, col = "white", nlevels = 8)
    if (best$p == pv)
      points(best$w_min, best$w_max, pch = 4, col = "red", cex = 2, lwd = 3)
  }
  
  invisible(list(best = best, grid = grid))
}

if(exists("validation")) {
alt_fit <- fit_alt_weights(validation)
alt_fit <- fit_alt_weights(
  validation,
  min_grid = seq(0.1, 0.7, by = 0.05),   # below the static 0.7
  max_grid = seq(0.7, 0.95, by = 0.05),   # at/above the static 0.7
  p_grid   = c(1, 2, 3)
)
}


fit_combined_weights <- function(validation,
                                 min_grid = seq(0.20, 0.60, by = 0.05),
                                 max_grid = seq(0.75, 0.95, by = 0.05),
                                 p_grid   = c(1, 2, 3),
                                 tau_grid = c(0.5, 1.0, 1.5, 2.0)) {
  if (!USE_AMBIENT) { warning("USE_AMBIENT is FALSE"); return(NULL) }
  
  fitdat <- validation %>%
    filter(!is.na(shadow_fraction), !is.na(ambient_fraction),
           !is.na(light_mol), !is.na(surface_par),
           !is.na(water_depth_m), !is.na(sun_altitude), !is.na(Kd))
  
  message("Tuning combined alt+depth on ", nrow(fitdat), " rows")
  
  water_atten <- exp(-fitdat$Kd * fitdat$water_depth_m)
  sin_alt     <- sin(pmax(fitdat$sun_altitude, 0))
  tau         <- pmax(fitdat$Kd * fitdat$water_depth_m, 0)
  
  grid <- expand.grid(w_min = min_grid, w_max = max_grid,
                      p = p_grid, tau_s = tau_grid)
  grid <- grid[grid$w_max > grid$w_min, ]
  grid$rmse <- grid$r2_1to1 <- grid$bias <- NA_real_
  
  for (i in seq_len(nrow(grid))) {
    A <- (1 - sin_alt)^grid$p[i]
    D <- 1 - exp(-tau / grid$tau_s[i])
    combined <- 1 - (1 - A) * (1 - D)
    w_amb <- grid$w_min[i] + (grid$w_max[i] - grid$w_min[i]) * combined
    w_amb <- pmin(pmax(w_amb, 0), 1)
    eff   <- (1 - w_amb) * fitdat$shadow_fraction + w_amb * fitdat$ambient_fraction
    pred  <- fitdat$surface_par * eff * water_atten
    res   <- pred - fitdat$light_mol
    grid$rmse[i]    <- sqrt(mean(res^2, na.rm = TRUE))
    grid$bias[i]    <- mean(res, na.rm = TRUE)
    ss_res <- sum(res^2); ss_tot <- sum((fitdat$light_mol - mean(fitdat$light_mol))^2)
    grid$r2_1to1[i] <- 1 - ss_res/ss_tot
  }
  grid <- grid[order(grid$rmse), ]
  cat("\nBest combined (alt+depth) weights:\n")
  print(head(grid, 15), row.names = FALSE, digits = 3)
  invisible(grid)
}

if(exists("validation")){
fit_combined <- fit_combined_weights(validation)
}

# To apply the recommendation:
#   USE_ALT_DEPENDENT_WEIGHTS <- TRUE

# ── UTM PROJECTION ────────────────────────────────────────────────────────────
# all DEMs are reprojected to UTM zone 55S before processing
# EPSG:32755 = WGS84 UTM Zone 55S, correct for Lizard Island (~145°E, ~15°S)
# reprojecting puts x, y, and z all in metres, so zscale = 1 (no unit conversion needed)
# suncalc and tideharmonics still need WGS84 lat/lon — get_site_centroid() handles this
UTM_CRS <- "EPSG:32755"
#utm is better because it is exact, no approximation needed to get meter coordinates
#so all spatial operations are preformed in meters

# ── kd_method_label ───────────────────────────────────────────────────────────
# One-stop label describing WHICH Kd method produced this run's attenuation, so
# every plot title says it explicitly (satellite vs reference vs spectral) and
# satellite/reference figures can't be confused at a glance. Reads the same
# toggles the pipeline uses, so it always matches the run that made the figure.
kd_method_label <- function() {
  if (USE_CONSTANT_KD)    return(sprintf("constant Kd (%.2f)", KD_CONSTANT))
  if (USE_REFERENCE_KD)   return("reference Kd")
  if (USE_SPECTRAL_KD)    return("spectral Kd")
  base <- paste0("satellite Kd (", KDPAR_COL, ")")
  if (exists("KDPAR_OVERRIDE") && length(KDPAR_OVERRIDE) > 0)
    base <- paste0(base, " + overrides")
  base
}
KD_METHOD <- kd_method_label()   # evaluated once; re-call kd_method_label() if you flip toggles
print(KD_METHOD)

# ══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# reusable code chunk recipes that I made to be called later
#
# DEFINITION ORDER MATTERS: make_ambient is defined BEFORE prepare_dem because
# prepare_dem now calls it inline (so the ambient layer ships with the DEM
# bundle and only computes once, before the hourly loop).
# ══════════════════════════════════════════════════════════════════════════════

# ── make_ambient ──────────────────────────────────────────────────────────────
#  Compute the ambient occlusion matrix ONCE per DEM. Same orientation/units as
#  ray_shade so it samples through the exact same extract path.
#  zscale must match what you pass to ray_shade (= dem$zscale, UTM metres).
#  multicore = TRUE would parallelise across cores but needs the parallel
#  backend set up beforehand; leaving FALSE for portability.
#
#  CHANGE: maxsearch is now derived per-DEM from AMBIENT_MAXSEARCH_M (a
#  physical distance in METRES). pixel_size on this DEM is `zscale`, so the
#  conversion is simply maxsearch_m / zscale. This keeps the physical ray
#  distance the same on every site even when native pixel sizes differ.
make_ambient <- function(dem_matrix, zscale,
                         sunbreaks   = AMBIENT_SUNBREAKS,
                         maxsearch_m = AMBIENT_MAXSEARCH_M) {
  maxsearch_px <- max(1, round(maxsearch_m / zscale))
  message("  Computing ambient occlusion (sunbreaks=", sunbreaks,
          ", maxsearch=", maxsearch_px, " px = ", maxsearch_m, " m) ...")
  amb <- ambient_shade(
    heightmap = dem_matrix,
    sunbreaks = sunbreaks,
    maxsearch = maxsearch_px,
    zscale    = zscale,
    multicore = FALSE
  )
  message("    ambient range: ", round(min(amb, na.rm = TRUE), 3),
          " to ",               round(max(amb, na.rm = TRUE), 3),
          " | mean: ",          round(mean(amb, na.rm = TRUE), 3))
  amb
}

# ── prepare_dem ───────────────────────────────────────────────────────────────
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
#      above) and converted to pixels per-DEM, so ray distances are physically
#      consistent across sites.
prepare_dem <- function(site, year,
                        reef_radius_px = REEF_RADIUS_PX_TARGET) {
  path <- file.path(DEM_DIR, paste0(site, year, ".tif"))
  message("  Loading DEM: ", path)
  dem  <- rast(path) |> trim()
  message("  Raw DEM CRS: ", crs(dem, describe=TRUE)$name)
  
  # reproject to UTM zone 55S so x, y, z are all in metres
  # this means zscale is meaningful (1 horizontal metre = 1 vertical metre)
  dem  <- project(dem, UTM_CRS)
  message("  Reprojected to UTM 55S (EPSG:32755)")
  message("  DEM extent (metres): ", paste(round(as.vector(ext(dem))), collapse = " "))
  
  # ── CROP TO REEF DISC BOUNDING BOX (before aggregate, for speed) ──────────
  # The raw DEM extent can be much wider than the 6.5 m reef (e.g. CB19 covers
  # ~23 m square). After deep-fill those off-reef pixels contribute nothing,
  # but rayshader still has to march rays across them. Cropping to the disc
  # bbox + a small margin shrinks the matrix without changing physics.
  CROP_MARGIN_M <- 0.5   # cosmetic buffer around the 6.5 m disc
  
  # data centroid on the native raster (mean x/y of non-NA cells)
  xy_native <- terra::xyFromCell(dem, which(!is.na(values(dem))))
  if (nrow(xy_native) == 0)
    stop("prepare_dem: DEM has no non-NA pixels after reproject for ", site, year)
  cx <- mean(xy_native[, 1]); cy <- mean(xy_native[, 2])
  
  half <- REEF_RADIUS_M + CROP_MARGIN_M
  crop_ext <- terra::ext(cx - half, cx + half, cy - half, cy + half)
  crop_ext <- terra::intersect(crop_ext, terra::ext(dem))   # stay inside data
  dem <- terra::crop(dem, crop_ext)
  message("  Cropped to reef bbox: ",
          paste(round(as.vector(ext(dem)), 2), collapse = " "),
          " (", round(ncol(dem)), " x ", round(nrow(dem)), " px)")
  
  # ── native pixel size on the reprojected (un-aggregated) DEM ─────────────
  res_x          <- (ext(dem)$xmax - ext(dem)$xmin) / ncol(dem)
  res_y          <- (ext(dem)$ymax - ext(dem)$ymin) / nrow(dem)
  native_px_size <- (res_x + res_y) / 2          # metres per native pixel
  
  # ── per-DEM downsample factor — targets REEF DISC resolution, not matrix ──
  # We want the 6.5 m reef to occupy ~reef_radius_px pixels on the downsampled
  # grid, so the disc has the same pixel density on every site. Solve for the
  # aggregation factor:
  #     target_px_size = REEF_RADIUS_M / reef_radius_px       (m per downsampled px)
  #     fact           = target_px_size / native_px_size
  # Floor at MIN_DOWNSAMPLE so a DEM already at higher resolution than the
  # target keeps its native pixels (we never upsample).
  target_px_size <- REEF_RADIUS_M / reef_radius_px
  fact           <- max(MIN_DOWNSAMPLE, round(target_px_size / native_px_size))
  message("  Native size: ", nrow(dem), " x ", ncol(dem),
          " | native px=", round(native_px_size, 5), " m",
          " | target px=", round(target_px_size, 5), " m",
          " -> downsample factor ", fact)
  
  dem_small  <- aggregate(dem, fact = fact)
  
  # zscale = downsampled pixel size in metres (since x, y, z are all UTM metres)
  pixel_size <- native_px_size                    # convenience alias for the messages below
  zscale     <- native_px_size * fact
  message("  Pixel size (UTM metres): ", round(pixel_size, 6))
  message("  Zscale (downsampled): ",    round(zscale, 6),
          "  -> reef disc ~", round(REEF_RADIUS_M / zscale), " px radius")
  
  # ── ALL masking/filling done on the RASTER, so dem_small and the matrix
  #    derived from it share one identical footprint, extent and orientation.
  #
  # ──────────────────────────────────────────────────────────────────────────
  #  REEF MASK — fixed physical radius in METRES (replaces the 95th-percentile
  #  data-driven radius that was here before).
  #
  #  WHY FIXED, NOT DATA-DRIVEN
  #  --------------------------
  #  The previous adaptive radius was sized to "where the DEM happens to have
  #  non-NA pixels", which varies by photogrammetry batch (year 17 vs 18 vs
  #  19). Wider DEMs (CB19, RS19) got bigger discs; tightly-clipped DEMs
  #  (TB18, CB18) got smaller ones. Different physical reef areas across
  #  sites = unfair cross-site comparison. Switching to a fixed REEF_RADIUS_M
  #  gives every site the SAME 13 m disc, regardless of how the DEM was clipped.
  #
  #  WHAT THIS DOES
  #  --------------
  #  1. Find where the DEM ACTUALLY has data (data_present = !is.na).
  #  2. Centre the disc on the DATA centroid, not the bbox centre — robust to
  #     off-centre surveys (GT19 etc.).
  #  3. Convert REEF_RADIUS_M (metres) to a pixel radius for THIS DEM. zscale
  #     IS the pixel size in m on dem_small, so rad_px = REEF_RADIUS_M / zscale.
  #  4. Intersect the disc with data_present. KEY LINE: a pixel is "inside"
  #     only if it is BOTH within the disc AND a real data pixel. This means:
  #       • Tightly clipped DEM smaller than 13 m: disc clips to the available
  #         data (no fake pixels invented).
  #       • Unclipped DEM with rotated-square fill outside the reef: disc
  #         cleanly trims that off — fill cannot leak into the reef region.
  #       • Rectangular DEM (GT19): disc clips to the rectangular footprint,
  #         which is still consistent across sites because the data centroid
  #         is what we centre on.
  # ──────────────────────────────────────────────────────────────────────────
  
  nr  <- nrow(dem_small)
  nc  <- ncol(dem_small)
  
  # values() is ROW-MAJOR; build a row-major logical of where real data sits,
  # then reshape to an [nr x nc] matrix matching row()/col() indexing.
  v0           <- values(dem_small)
  data_present <- matrix(!is.na(v0), nrow = nr, ncol = nc, byrow = TRUE)
  
  if (!any(data_present))
    stop("prepare_dem: DEM has no non-NA pixels after reproject/trim for ",
         site, year)
  
  rr <- row(data_present)
  cc <- col(data_present)
  
  # data centroid (mean row/col of real pixels) — robust to off-centre surveys
  row_c <- mean(rr[data_present])
  col_c <- mean(cc[data_present])
  
  # convert physical radius to pixel radius for this DEM's downsampled grid
  rad_px <- REEF_RADIUS_M / zscale
  
  dist_all    <- sqrt((rr - row_c)^2 + (cc - col_c)^2)
  circle_mask <- dist_all <= rad_px
  
  # KEY INTERSECTION: inside == within the disc AND genuinely real data.
  inside_mat  <- circle_mask & data_present
  
  if (!any(inside_mat))
    stop("prepare_dem: fixed-metres disc + data_present intersection is empty ",
         "for ", site, year, ". Check REEF_RADIUS_M and DEM extent.")
  
  message("  Reef mask: ", sum(inside_mat), " px inside ",
          "(centroid row/col ", round(row_c), "/", round(col_c),
          ", radius ", round(rad_px), " px = ", REEF_RADIUS_M, " m)")
  
  # write back: everything outside the intersected mask -> NA (row-major vector)
  v <- v0
  v[!as.vector(t(inside_mat))] <- NA
  values(dem_small) <- v
  
  # ──────────────────────────────────────────────────────────────────────────
  #  DEEP-FILL BEFORE RAYTRACING (replaces the old reef_min − relief fill)
  #
  #  WHY 100 M BELOW REEF_MIN
  #  ------------------------
  #  ambient_shade fires a fan of rays from each pixel and tests what fraction
  #  of the sky those rays can reach without hitting terrain. If off-reef
  #  pixels are at heights comparable to reef_min (which they could be with
  #  the old fill = reef_min − relief on a low-relief site, or where a tilted
  #  rotated-square fill happens to interpolate to a high corner), the rays
  #  would see them as terrain and report phantom occlusion at the disc edge.
  #
  #  Setting the off-disc fill to reef_min − 100 m guarantees those pixels are
  #  far enough below the reef that no ambient ray from inside the disc can
  #  reach them within a sensible vertical search. The off-reef region is now
  #  invisible to both ambient_shade and ray_shade — the reef interior is
  #  computed as if the disc is the entire DEM, just embedded in a deep ocean
  #  floor of "empty" cells.
  #
  #  Crucially, this happens BEFORE make_ambient() is called below.
  # ──────────────────────────────────────────────────────────────────────────
  reef_min <- min(values(dem_small), na.rm = TRUE)
  fill_val <- reef_min - 100   # 100 m below; well outside any sensible search
  dem_small[is.na(dem_small)] <- fill_val
  
  # matrix derived LAST, from the filled raster (rayshader needs real numbers everywhere)
  dem_matrix <- raster_to_matrix(raster(dem_small))
  
  # keep the reef mask so we can blank the border to NA AFTER raytracing.
  # inside_mat is already [nr x nc] in row/col indexing; transpose to match the
  # rayshader/matrix orientation (same convention as the shadow matrices).
  reef_mask_mat <- t(inside_mat)
  
  message("  DEM loaded: ", nrow(dem_matrix), " x ", ncol(dem_matrix),
          " | raster: ", nrow(dem_small), " x ", ncol(dem_small))
  
  # ── AMBIENT OCCLUSION — computed ONCE here (geometry-only, sun-independent) ──
  # Same matrix orientation/units as ray_shade so it samples through the exact
  # same extract path. zscale matches what we pass ray_shade.
  amb_matrix <- make_ambient(dem_matrix, zscale)
  
  list(matrix = dem_matrix, rast = dem_small, zscale = zscale,
       reef_mask = reef_mask_mat,
       ambient   = amb_matrix)         # diffuse occlusion layer
}
#returns: a matrix for rayshader, a raster for coordinate extraction,
#the zscale, the reef mask, and the ambient occlusion matrix.

# ── load_hobo_locations ───────────────────────────────────────────────────────
# build a sf point object from geom_new coordinates for a given site/year.
# coords are in the local transverse mercator projection used when the
# geom script extracted them from the DEMs. We need to know that projection
# explicitly so st_transform can convert to UTM later.
load_hobo_locations <- function(site, year, geom) {
  tm_crs <- paste0("+proj=tmerc +lat_0=-14.6989736557007", #this defines the coordinate system the numbers are in
                   " +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0",
                   " +datum=WGS84 +units=m +no_defs")
  locs <- geom %>%
    filter(site == !!site, year == !!year) %>%
    filter(!is.na(lon), !is.na(lat)) %>%
    distinct(hoboID, .keep_all = TRUE)
  
  if (nrow(locs) == 0) {
    stop("No logger locations found in geom for site=", site, " year=", year)
  }
  
  message("  Loaded ", nrow(locs), " logger locations from geom for ",
          site, year)
  st_as_sf(locs, coords = c("lon", "lat"), crs = tm_crs)
}

# ── get_sun ───────────────────────────────────────────────────────────────────
#this uses a timestamp, lat and long to get the sun altitude (height above horizon in radians)
#and azimuth, which is the compass direction in radians
#suncalc requires WGS84 lat/lon — get_site_centroid() reprojects from UTM for this
get_sun <- function(time, lat, lon) {
  sun <- getSunlightPosition(date = time, lat = lat, lon = lon)
  list(altitude = sun$altitude, azimuth = sun$azimuth)
}
#this determines the angle / strength of the shadows

# ── make_shadow ───────────────────────────────────────────────────────────────
#this runs the raytracing in the background, sending the dem_matrix that I made with prepare_dem
#and sun positioning to do raytracing using rayshader, this returns a grid of values
#0 is fully shaded, and 1 is fully lit. If the sun is below the horizon (or
#very near it) it skips the measurements to save time AND because near-horizon
#rays produce unreliable shadows on noisy DEMs (long shadows triggered by
#tiny coral bumps the photogrammetry happened to resolve).
#
#  CHANGE: maxsearch is now derived per-DEM from MAXSEARCH_M (a physical
#  distance in METRES). pixel_size on this DEM is `zscale`, so the conversion
#  is simply MAXSEARCH_M / zscale. This keeps direct-beam ray distances the
#  same on every site regardless of native pixel size. At TB18's zscale
#  (~0.0123 m/px on a 1062px DEM) this is ~700 px, matching the old
#  MAXSEARCH_FACT=700 — so TB18 reproduces previous shadows and other sites
#  become consistent with it.
make_shadow <- function(dem_matrix, sun, zscale) {
  if (sun$altitude <= 0.05) {
    message("    (sun below horizon — skipping raytrace)")
    return(NULL)
  }
  
  maxsearch_px <- max(1, round(MAXSEARCH_M / zscale))
  
  shadow <- ray_shade(
    dem_matrix,
    sunaltitude = sun$altitude * 180 / pi, #gets sun altitude in degrees from radians
    sunangle    = sun$azimuth  * 180 / pi, #gets sun azimuth in degrees from radians
    zscale      = zscale,
    maxsearch   = maxsearch_px  # MAXSEARCH_M metres, converted to px via zscale
  )
  
  message("    shadow range: ", round(min(shadow, na.rm = TRUE), 3),
          " to ",               round(max(shadow, na.rm = TRUE), 3),
          " | mean: ",          round(mean(shadow, na.rm = TRUE), 3),
          " | maxsearch=",      maxsearch_px, " px (", MAXSEARCH_M, " m)")
  shadow
}


# ── calc_water_attenuation ────────────────────────────────────────────────────
# Beer-Lambert law: exp(−Kd × depth). Single-KdPAR path (default).
# Returns the fraction of surface PAR surviving to depth; deeper logger = smaller.
# Guard: if depth <= 0 (sensor briefly above water at low tide) or Kd is NA,
# return NA rather than producing a spurious value. The downstream pred_par will
# be NA for that hour/sensor, which the validation join filters out.
calc_water_attenuation <- function(Kd, water_depth_m) {
  ifelse(
    is.na(Kd) | water_depth_m <= 0,
    NA_real_,
    exp(-Kd * water_depth_m)
  )
}

# ── calc_water_attenuation_spectral ───────────────────────────────────────────
# Depth-dependent spectral path. Instead of one KdPAR, attenuate each band at its
# own Kd and sum the surviving PAR fraction at depth; depth-dependence emerges
# automatically (Canto method). f_lambda = surface photon weight per band (∝ λ,
# same assumption as kdpar_photon; Kirk 2011). Drops NA bands and renormalises so
# a site missing one band (e.g. CB18, 9 bands) still computes.
calc_water_attenuation_spectral <- function(kd_bands, waves, water_depth_m) {
  if (is.na(water_depth_m) || water_depth_m <= 0) return(NA_real_)   # depth guard
  ok <- is.finite(kd_bands)
  if (sum(ok) < 2) return(NA_real_)   # need at least 2 bands for a meaningful sum
  kb <- kd_bands[ok]; wv <- waves[ok]
  f  <- wv / sum(wv)
  sum(f * exp(-kb * water_depth_m))
}


# ── derive_reference_kd ───────────────────────────────────────────────────────
# Takes a VECTOR of reference logger IDs. For each, inverts Beer-Lambert per
# hour; then aggregates across refs via REF_AGG_FUN (median by default).
# Returns a named vector: hourly POSIXct (as character) -> aggregated Kd.
#
# Per-logger inversion:
#   Kd(t) = -ln( ref_light(t) / surface_par(t) ) / ref_depth(t)
# Guards: needs ref_light < surface_par (transmittance in (0,1) — otherwise
# log() is undefined or positive) AND depth > 0 (otherwise division blows up).
#
# WHY MULTIPLE REFS: a single reference logger that happens to be slightly
# shaded at some hour produces a transmit value that's too low, which makes
# -log(transmit) too high, which spikes Kd at that hour, which contaminates
# every other logger's prediction at that hour. Median across 3 refs is robust
# to one ref blowing up — the middle value is reported regardless. Won't help
# if all 3 share the same bias (see the finder docs below).
#
# WHY MEDIAN NOT MEAN: median ignores one bad ref. Mean would drag the
# aggregate toward the bad ref. Switch via REF_AGG_FUN if you want to
# experiment.
derive_reference_kd <- function(env_site, ref_ids, cap_tide_m,
                                agg_fun = REF_AGG_FUN) {
  if (length(ref_ids) == 0) return(NULL)
  
  # per-ref: list of named vectors (timestamp -> Kd_t), one element per ref.
  # Same shape as the old single-ref output, just stored side-by-side rather
  # than overwriting.
  per_ref <- lapply(ref_ids, function(ref_id) {
    ref <- env_site %>%
      filter(as.character(hoboID) == ref_id) %>%
      mutate(datetime_pos = floor_date(datetime_pos, TIMESTEP_INTERVAL)) %>%
      group_by(datetime_pos) %>%
      summarise(surface_par = mean(surface_par, na.rm = TRUE),
                ref_light   = mean(light_mol,   na.rm = TRUE),
                tide_cm     = mean(tide_cm,     na.rm = TRUE),
                dpHB        = mean(dpHB,        na.rm = TRUE), .groups = "drop") %>%
      mutate(ref_depth = (tide_cm / 100) - dpHB - cap_tide_m,
             transmit  = ref_light / surface_par,
             Kd_t      = ifelse(transmit > 0 & transmit < 1 & ref_depth > 0,
                                -log(transmit) / ref_depth, NA_real_))
    
    # diagnostic — eyeball this. transmit > 1 or wildly variable Kd means the
    # ref logger isn't actually open / has a unit mismatch with surface_par /
    # has a wrong depth — DON'T use it without investigating.
    message("  [ref-Kd] ", ref_id,
            " | transmit ",
            round(min(ref$transmit, na.rm = TRUE), 3), "-",
            round(max(ref$transmit, na.rm = TRUE), 3),
            " | Kd(t) ",
            round(min(ref$Kd_t, na.rm = TRUE), 3), "-",
            round(max(ref$Kd_t, na.rm = TRUE), 3),
            " | median ", round(median(ref$Kd_t, na.rm = TRUE), 3))
    
    setNames(ref$Kd_t, as.character(ref$datetime_pos))
  })
  
  # union of timestamps across all refs, then aggregate per timestamp.
  # If ref 1 covers hours 11-17 but ref 2 covers 12-18, all_ts is 11-18; hours
  # where a ref has no value come back as NA via v[all_ts] and are ignored by
  # the na.rm = TRUE in REF_AGG_FUN.
  all_ts <- sort(unique(unlist(lapply(per_ref, names))))
  Kd_mat <- vapply(per_ref,
                   function(v) v[all_ts],
                   numeric(length(all_ts)))
  if (!is.matrix(Kd_mat)) Kd_mat <- matrix(Kd_mat, ncol = length(per_ref))
  Kd_agg <- apply(Kd_mat, 1, agg_fun)
  names(Kd_agg) <- all_ts
  
  # how many refs actually contributed each hour — if this is way less than
  # length(ref_ids) on average, your refs are frequently hitting the
  # transmit/depth guards (i.e. they aren't all good refs).
  n_contrib <- apply(Kd_mat, 1, function(x) sum(is.finite(x)))
  message("  [ref-Kd] aggregated across ", length(ref_ids), " refs (",
          paste(ref_ids, collapse = ", "), ")",
          " | median Kd ", round(median(Kd_agg, na.rm = TRUE), 3),
          " | range ", round(min(Kd_agg, na.rm = TRUE), 3), "-",
          round(max(Kd_agg, na.rm = TRUE), 3),
          " | mean refs/hour: ", round(mean(n_contrib), 2),
          " of ", length(ref_ids))
  
  Kd_agg
}



# ── extract_hobo_dem ──────────────────────────────────────────────────────────
#this looks up hobo logger locations and stores it as a new column called dem elevation, which is not used in the final prediction but good to have
#reproject hobo points to UTM to match the DEM (now in UTM after reprojection in prepare_dem)
extract_hobo_dem <- function(dem_rast, hobo_sf) {
  hobo_utm         <- st_transform(hobo_sf, UTM_CRS)
  vals             <- terra::extract(dem_rast, vect(hobo_utm))
  hobo_sf$dem_elev <- vals[, 2]
  hobo_sf
}

# ── extract_hobo_shadow ───────────────────────────────────────────────────────
#this takes the shadow grid from make_shadow, and puts it in to a spatial raster with the same coordinates as the raster, than samples each loggers exact location
#to get the logger's shadow fraction for that hour
#the shadow raster is now in UTM (because the DEM is in UTM), so hobo points are reprojected to UTM
extract_hobo_shadow <- function(shadow_matrix, dem_rast, hobo_sf) {
  shadow_rast <- rast(t(shadow_matrix))
  ext(shadow_rast) <- ext(dem_rast)
  crs(shadow_rast) <- crs(dem_rast)   # UTM 55S
  
  # reproject hobo points to UTM to match the shadow raster
  # UTM units are metres so buffer distance is directly in metres — no degree conversion needed
  hobo_utm <- st_transform(hobo_sf, UTM_CRS)
  
  # buffer in meters: BUFFER_M = 0.01 = 1 cm. Small buffer averages over a
  # tiny disc around the logger point so we don't get fooled by a single
  # noisy pixel. If logger placement uncertainty is larger (e.g. 5 cm), bump
  # BUFFER_M to match.
  hobo_buffered <- st_buffer(hobo_utm, dist = BUFFER_M)
  
  # extract mean shadow value within each buffer
  vals <- terra::extract(shadow_rast, vect(hobo_buffered), fun = mean, na.rm = TRUE)
  vals[, 2]
}

# ── extract_hobo_ambient ──────────────────────────────────────────────────────
#  Sample the ambient occlusion matrix at each logger, using identical machinery
#  to extract_hobo_shadow so the two layers line up pixel-for-pixel and any
#  buffer/extract changes apply to both.
extract_hobo_ambient <- function(ambient_matrix, dem_rast, hobo_sf) {
  amb_rast <- rast(t(ambient_matrix))
  ext(amb_rast) <- ext(dem_rast)
  crs(amb_rast) <- crs(dem_rast)            # UTM 55S
  
  hobo_utm      <- st_transform(hobo_sf, UTM_CRS)
  hobo_buffered <- st_buffer(hobo_utm, dist = BUFFER_M)   # match extract_hobo_shadow buffer
  
  vals <- terra::extract(amb_rast, vect(hobo_buffered), fun = mean, na.rm = TRUE)
  vals[, 2]
}

# ── get_site_centroid ─────────────────────────────────────────────────────────
# finds the center point of the DEM in WGS84 lat/lon degrees (EPSG:4326)
# even though the DEM is now in UTM, this reprojects back to WGS84 for suncalc
# suncalc requires geographic lat/lon — UTM coordinates are meaningless to it
get_site_centroid <- function(dem_rast) {
  e <- ext(project(dem_rast, "EPSG:4326"))
  list(lat = (e$ymin + e$ymax) / 2, lon = (e$xmin + e$xmax) / 2)
}


# ══════════════════════════════════════════════════════════════════════════════
#  PER-TIMESTEP PROCESSING
#the core calculations for each timestep (1 hour by default)
#  Blends a directional (ray_shade) and a diffuse (ambient_shade) layer with
#  weights that can be static or altitude-dependent.
# ══════════════════════════════════════════════════════════════════════════════
process_timestep <- function(env_block, dem_matrix, dem_rast, dem_zscale,
                             hobo_sf, centroid, dem_reef_mask, dem_ambient,
                             kd_bands = NULL) {
  
  # 1. calculate solar position for the chosen hour at the chosen site
  # centroid is in WGS84 lat/lon (reprojected from UTM by get_site_centroid)
  timestamp <- env_block$datetime_pos[1]
  sun <- get_sun(timestamp, centroid$lat, centroid$lon)
  
  message("  >> Raytracing: ", format(timestamp, "%Y-%m-%d %H:%M"),
          " | sun alt: ", round(sun$altitude * 180 / pi, 1), "\u00b0",
          " | sun az: ",  round(sun$azimuth  * 180 / pi, 1), "\u00b0",
          " | Kd: ",      round(env_block$Kd[1], 4))
  
  # 2. directional shadow (sun beam) from the dem
  # zscale already absorbs the UTM-metres units
  shadow <- make_shadow(dem_matrix, sun, zscale = dem_zscale)
  # NA-mask the off-reef border for sampling/plotting. With the deep-fill in
  # prepare_dem, the off-reef pixels were already 100 m below the reef during
  # raytracing so they could not have affected reef-interior shadows — this
  # mask is now mainly cosmetic (clean disc edges in plots, no out-of-disc
  # extract values). (guard: shadow is NULL at night)
  if (!is.null(shadow)) shadow[!dem_reef_mask] <- NA
  
  # 2b. ambient occlusion (diffuse) — sun-independent, blank border the same way
  # so the two layers have an identical sampling footprint
  ambient_masked <- dem_ambient
  ambient_masked[!dem_reef_mask] <- NA
  
  # 3. sample BOTH layers at each hobo logger.
  # if the sun is below the horizon (dark hours) shadow is NA and skipped;
  # ambient is geometry-only so it is always sampled, but it doesn't matter
  # because surface_par is ~0 in dark hours anyway -> pred_par ~0.
  shadow_vals <- if (!is.null(shadow)) {
    extract_hobo_shadow(shadow, dem_rast, hobo_sf)
  } else {
    rep(NA_real_, nrow(hobo_sf))
  }
  ambient_vals <- extract_hobo_ambient(ambient_masked, dem_rast, hobo_sf)
  
  message("    shadow per logger:  ",
          paste(hobo_sf$hoboID, "=", round(shadow_vals, 3), collapse = ", "))
  message("    ambient per logger: ",
          paste(hobo_sf$hoboID, "=", round(ambient_vals, 3), collapse = ", "))
  
  # 4. weights for THIS hour.
  #    USE_AMBIENT=FALSE: pure direct-beam (w_amb=0, eff_light=shadow).
  #    USE_AMBIENT=TRUE + USE_ALT_DEPENDENT_WEIGHTS=TRUE: physically motivated
  #         altitude-dependent W_AMBIENT (rises at low sun).
  #    USE_AMBIENT=TRUE + USE_ALT_DEPENDENT_WEIGHTS=FALSE: static W_AMBIENT.
  #OLD: must now vary per logger
  #w_amb_now <- if (USE_AMBIENT) {
  #  if (USE_ALT_DEPENDENT_WEIGHTS) w_ambient_alt(sun$altitude) else W_AMBIENT
  #} else 0
  # NEW: defer per-row computation into the mutate below
  sun_alt_rad <- sun$altitude   # scalar, same for all loggers this hour
  
  message("    sun_alt=", round(sun_alt_rad*180/pi, 1), "°  ",
          "(w_amb varies per logger by depth)")
  
  # 5. pair fractions with logger metadata (this hour, this site),
  #    join env data (depth, surface par, Kd), blend the directional + diffuse
  #    layers, and compute the final predicted value as
  #      surface_par × (w_direct·shadow + w_ambient·ambient) × water_atten
  shadow_df <- data.frame(
    hoboID_env       = hobo_sf$hoboID,
    shadow_fraction  = shadow_vals,
    ambient_fraction = ambient_vals
  )
  
  result <- shadow_df %>%
    left_join(
      env_block %>%
        mutate(hoboID_env = as.character(hoboID)) %>%
        select(hoboID_env, Kd, water_depth_m, surface_par, site, year),
      by = "hoboID_env"
    ) %>%
    mutate(
      datetime_pos = timestamp,
      sun_altitude = sun_alt_rad,
      sun_azimuth  = sun$azimuth,
      
      # per-logger weights — depend on this row's depth + Kd
      w_ambient = if (USE_AMBIENT) {
        w_ambient_combined(sun_alt_rad, water_depth_m, Kd)
      } else 0,
      w_direct = 1 - w_ambient,
      
      water_atten = if (USE_SPECTRAL_KD && !is.null(kd_bands)) {
        vapply(water_depth_m,
               function(z) calc_water_attenuation_spectral(kd_bands, WAVES_PAR, z),
               numeric(1))
      } else {
        calc_water_attenuation(Kd, water_depth_m)
      },
      eff_light = if (USE_AMBIENT) {
        w_direct * shadow_fraction + w_ambient * ambient_fraction
      } else {
        shadow_fraction
      },
      pred_par = surface_par * eff_light * water_atten
    )
  
  as_tibble(result)
}


# ══════════════════════════════════════════════════════════════════════════════
#  main code loop that puts it all together
# master function that does following:
#
# loops over every site x year combination in data
# loads elevation map and logger positions ONCE
# builds a schedule of hourly daytime timesteps to process
# and for every hour that qualifies (sun > 20° over horizon), raytraces DEM and
# computes PAR at every logger location
#
# the outer loop sets up the machinery for every site/year, and the inner loop
# runs it one hour at a time
#
# inputs:  env_long  - full environmental data table
#          geom      - logger coordinate lookup table
#          dem_store - (optional) pre-loaded DEMs/hobo points to skip prepare_dem
#
# output: one table with one row per logger per timestep, with pred_par
# ══════════════════════════════════════════════════════════════════════════════

# ── helper: one canonical site/year key, used everywhere ──────────────────────
# Defined BEFORE run_pipeline because run_pipeline calls it. Keys look like
# "TB18", "CB19" — match the format in kdpar_lookup, capture_tide, and
# REFERENCE_LOGGER so joins/lookups all line up.
make_key <- function(site, year) paste0(site, year)   # "TB" + 18 -> "TB18"

run_pipeline <- function(env_long, geom, dem_store = NULL) {
  #results starts as an empty list and fills one chunk at a time.
  #dem_store (optional): named list keyed "TB18" etc, built ONCE above. If
  #supplied, the spatial data (DEM, ambient occlusion, hobo points) is REUSED
  #from it instead of calling prepare_dem() again — prepare_dem (esp.
  #ambient_shade) is the slow step and was previously run twice per site/year.
  
  results <- list()
  #get the unique site codes present (CB, L1, GT, etc)
  sites   <- unique(env_long$site)
  #unique skips duplicates so we only do one site once
  
  message("=================================================")
  message("Starting pipeline")
  message("Sites to process: ", paste(sites, collapse = ", "))
  message("=================================================")
  
  #outer loop 1 -----------------------------------------------------------------
  #each site has its own dem file and logger positions, so it processes them one at a time
  for (s in sites) {
    years <- unique(env_long$year[env_long$site == s])
    #find which years of data exist for each site
    #the env_long$year[..... ] means give me the year column but for only rows where the site column equals s
    
    #outer loop 2 -----------------------------------------------------------------
    #---------step1 load spatial data
    #iterate over years within site:
    #load spatial data once per site/year
    #when they are outside hourly loop we only do this once, since it takes a lot of time (it can hang R if done inside the loop)
    for (y in years) {
      message("\n--- Site: ", s, " | Year: ", y, " ---")
      key <- make_key(s, y)
      #load dem file, downsample it, convert to matrix
      #returns $matrix, $rast, $zscale, $reef_mask, $ambient
      # ── Reuse cached spatial data (from dem_store) if available ────────────
      # The dem_store loop above already ran prepare_dem() + extract_hobo_dem()
      # once for every site/year in the run. Reuse those instead of rebuilding
      # them here (prepare_dem, esp. ambient_shade, is the slow step). Falls
      # back to loading from scratch if a key is missing or dem_store is NULL.
      if (!is.null(dem_store) && !is.null(dem_store[[key]])) {
        message("  Reusing cached DEM/hobo for ", key, " (skipping prepare_dem)")
        dem      <- dem_store[[key]]$dem      # $matrix, $rast, $zscale, $reef_mask, $ambient
        hobo_sf  <- dem_store[[key]]$hobo_sf  # already has extract_hobo_dem() applied
      } else {
        message("  No cached DEM for ", key, " — loading from scratch")
        dem      <- prepare_dem(s, y) #DEM is reprojected to UTM inside prepare_dem; ambient computed here too
        hobo_sf  <- load_hobo_locations(s, y, geom)   # uses geom_new, not shapefiles
        hobo_sf  <- extract_hobo_dem(dem$rast, hobo_sf) #sample DEM elevation at logger locations
      }
      centroid <- get_site_centroid(dem$rast) #reprojects UTM dem back to WGS84 lat/lon for suncalc
      message("  Site centroid: lat ", round(centroid$lat, 4),
              " lon ", round(centroid$lon, 4))
      
      # ── capture-tide for THIS site/year (replaces an old hardcoded "- 2") ──
      # one value, looked up from the capture_tide table built near the top.
      # If a site/year is missing from capture_tide (e.g. you added it without
      # updating the photogrammetry table), warn and fall back to -2 m so the
      # run doesn't silently NA out every depth.
      cap_tide_m <- capture_tide$capture_tide_m[
        capture_tide$site == s & capture_tide$year == y
      ]
      if (length(cap_tide_m) == 0 || is.na(cap_tide_m[1])) {
        warning("No capture_tide_m for site=", s, " year=", y,
                " — falling back to the old -2 m constant.")
        cap_tide_m <- 2
      } else {
        cap_tide_m <- cap_tide_m[1]
      }
      message("  Capture-tide correction: -", round(cap_tide_m, 3), " m",
              " (was a hardcoded -2 m)")
      
      # ── reference-logger time-varying Kd for THIS site (if toggled) ──────────
      # ref_kd_lookup is a named vector hour->Kd; the hourly_times mutate below
      # looks up each hour and uses the ref-derived Kd if available, else falls
      # back to the satellite KdPAR for the site.
      ref_kd_lookup <- NULL
      if (USE_REFERENCE_KD) {
        ref_ids <- REFERENCE_LOGGER[[make_key(s, y)]]
        if (!is.null(ref_ids) && length(ref_ids) > 0) {
          env_site_full <- env_long %>% filter(site == s, year == y)
          ref_kd_lookup <- derive_reference_kd(env_site_full, ref_ids, cap_tide_m)
          message("  Using ", length(ref_ids), " reference logger(s) (",
                  paste(ref_ids, collapse = ", "), ") for time-varying Kd")
        } else {
          message("  No reference loggers for ", make_key(s, y),
                  " — falling back to satellite Kd")
        }
      }
      
      #---------step2 build hourly schedule of timesteps
      #the raw data is 5 minute intervals, but we only raytrace 1 hour to save computing.
      #the section here collapses 5 min data to one row per logger per hour, computes sun angle, and throws away hours where the sun is too low to mean anything.
      #%>% chains everything together in order.
      hourly_times <- env_long %>%
        
        #keep only rows belonging to this site year
        filter(site == s, year == y) %>%
        
        #sort chronologically for floor date to work
        arrange(datetime_pos) %>%
        
        #snap every 5 min timestamp to its hour
        #11:05, 11:10, ..., 11:55 all become 11:00, which makes rows share the same datetime
        mutate(datetime_pos = floor_date(datetime_pos, TIMESTEP_INTERVAL)) %>%
        
        #group rows by hour and logger id, so each hour/logger gets its own independent group
        group_by(datetime_pos, hoboID) %>%
        
        #collapse the 5-min readings to the HOURLY MEAN per logger, so the
        #predictor (surface_par, tide, depth inputs) is mean-of-hour — matching
        #the validation join, which averages light_mol over the same hour.
        group_by(datetime_pos, hoboID, site, year) %>%
        summarise(
          surface_par = mean(surface_par, na.rm = TRUE),
          tide_cm     = mean(tide_cm,     na.rm = TRUE),
          dpHB        = mean(dpHB,         na.rm = TRUE),
          .groups = "drop"
        ) %>%
        
        #assign each row its Kd:
        #  - default = satellite KdPAR for this site/year
        #  - if ref-Kd is on AND we have a value for this hour, override with ref Kd
        #  - if both are missing, fall back to KD_CONSTANT
        # The Kd from env_long is intentionally NOT used here — it's the
        # sensor-back-calculated one, which would be circular.
        mutate(
          Kd = {
            sat_kd <- KDPAR_BY_KEY[[paste0(s, sprintf("%02d", as.integer(y) %% 100))]]
            if (is.null(sat_kd) || is.na(sat_kd)) sat_kd <- KD_CONSTANT
            if (USE_REFERENCE_KD && !is.null(ref_kd_lookup)) {
              rk <- ref_kd_lookup[as.character(datetime_pos)]
              ifelse(is.na(rk), sat_kd, rk)   # fall back to satellite where ref is NA
            } else {
              sat_kd
            }
          },
          # tide-corrected water column above sensor (recomputed here, not upstream)
          # tide_cm/100 = tide above datum (m); dpHB = sensor bed elevation (m, negative down)
          # cap_tide_m  = tide at this site's photogrammetry capture time,
          #               predicted by harmonic fit (replaces the old flat -2).
          water_depth_m = (tide_cm / 100) - dpHB - cap_tide_m
        ) %>%
        
        # Each site/year now carries its own KdPAR. The message lines that print
        # KD_CONSTANT will still print the fallback value, which is fine — the actual Kd
        # used is in env_block$Kd and prints via process_timestep's
        # "| Kd: round(env_block$Kd[1], 4)" message, so you can confirm per-site.
        
        #r applies the next mutate one row at a time; getSunlightPosition needs one single timestamp
        rowwise() %>%
        
        #compute solar altitude in degrees, store as sun_alt, converting radians -> degrees
        #centroid is in WGS84 lat/lon as required by suncalc
        mutate(
          sun_alt = getSunlightPosition(
            datetime_pos, centroid$lat, centroid$lon
          )$altitude * 180 / pi #return degrees
        ) %>%
        ungroup() %>%
        
        # discard hours where the sun is less than 20 degrees over the horizon.
        # at super low angles, the shadows hit DEM imperfections harder (long
        # shadows from tiny noise) AND the altitude-dependent W_AMBIENT is at
        # its highest, so the prediction relies heavily on the ambient layer.
        # we lose ~1 hour of early morning / late afternoon data per day.
        filter(sun_alt > 20) %>%
        
        #remove sun_alt column, since it served its purpose
        select(-sun_alt)
      
      # DIAGNOSTIC — confirm depth varies with tide; not strictly needed but
      # a useful eyeball check. swing = max - min depth over the deployment.
      print(
        hourly_times %>%
          group_by(hoboID) %>%
          summarise(
            depth_min  = round(min(water_depth_m), 2),
            depth_max  = round(max(water_depth_m), 2),
            depth_mean = round(mean(water_depth_m), 2),
            swing      = round(max(water_depth_m) - min(water_depth_m), 2)
          )
      )
      
      
      #get list of hours to raytrace, we only need one shadow map per hour
      timestamps <- unique(hourly_times$datetime_pos)
      
      #count how many steps we calculate, which is used for the message below
      n_steps    <- length(timestamps)
      
      message("  Timesteps to process: ", n_steps,
              " (", format(min(timestamps), "%Y-%m-%d %H:%M"),
              " to ", format(max(timestamps), "%Y-%m-%d %H:%M"), ")")
      
      #inner loop -----------------------------------------------------------------
      #---------step 3 inner loop, one raytracing per hour
      #for each hour: pull out all logger rows for that hour, pass them to
      #process_timestep, which raytraces once and extracts shadow + ambient values
      #for all loggers from that single shadow map.
      
      #seq_len(n_steps) produces 1,2,3,...,n_steps; safer than 1:n_steps because
      #if n_steps were 0, 1:0 would loop twice instead of zero times.
      for (i in seq_len(n_steps)) {
        
        #pick the i-th timestamp from the hours we decided
        t         <- timestamps[i]
        
        #pull out logger rows for this hour. process_timestep needs all loggers
        #at once so it can extract shadow values for each one from the shadow
        #map that was made earlier (one raytrace, many extractions).
        env_block <- hourly_times[hourly_times$datetime_pos == t, ]
        
        message("\n[", i, "/", n_steps, "] site: ", s, " | year: ", y,
                " | loggers in block: ", nrow(env_block),
                " | Kd: ", round(env_block$Kd[1], 4))
        
        #run the full prediction for the timestep
        #solar position -> raytrace + ambient sample -> Beer-Lambert -> predicted PAR
        # per-site 10-band Kd vector for the spectral path (NULL-safe if absent)
        site_kd_bands <- KD_BANDS_BY_KEY[[paste0(s, sprintf("%02d", as.integer(y) %% 100))]]
        
        res <- process_timestep(
          env_block,   #data frame: one row per logger for this hour
          dem$matrix,  #dem as a numeric matrix
          dem$rast,    #dem as a spatial object/raster for coordinate stuff (UTM)
          dem$zscale,  #zscale in UTM metres
          hobo_sf,     #logger locations as a spatial object
          centroid,    #site lat/lon (WGS84) for suncalc sun angle lookup
          dem$reef_mask,  # blank the off-reef border before sampling
          dem$ambient,    # precomputed ambient occlusion matrix
          site_kd_bands   # per-band Kd for spectral attenuation (NULL if not USE_SPECTRAL_KD)
        )
        
        #add hour's result tibble to results
        #length(results) + 1 puts it as the next empty slot
        results[[length(results) + 1]] <- res
      }
      
      message("\n  Done with site: ", s, " | year: ", y)
    }
  }
  
  message("\n=================================================")
  message("Pipeline complete :) Total result chunks: ", length(results))
  message("=================================================\n")
  
  #combine all hour's tibbles into one single df and return it
  #creating the final output of the data pipeline
  bind_rows(results)
}


# ══════════════════════════════════════════════════════════════════════════════
#  RUN — edit the `wanted` table; EVERYTHING below adapts to 1..N site/years
# ══════════════════════════════════════════════════════════════════════════════
#
#  HOW MULTI-SITE WORKS
#  ────────────────────
#  Nothing below hardcodes a single site. The flow is:
#    1. `wanted`     : the site/year combos you want this run.
#    2. dem_store    : DEMs/hobo points are loaded ONCE here and cached in a
#                      named list keyed "TB18", "CB19", ... so every plot reuses
#                      them instead of calling prepare_dem again.
#    3. run_pipeline : loops them internally, returns one pooled `output`.
#    4. per-site figs: PAR-per-logger panels and shadow maps LOOP over each
#                      site/year and save one file each.
#    5. pooled figs  : DLI + time-of-day scatter pool all rows; point SHAPE = the
#                      site/year, COLOUR = depth, and every sensor is labelled.
#
#  The one rule that makes it all work: a logger is only unique WITHIN a
#  site/year, so everywhere we used hoboID alone we now use a sitekey + hoboID.
# ══════════════════════════════════════════════════════════════════════════════

# ── FILTER DATA ───────────────────────────────────────────────────────────────
#  Edit ONLY this table. One row per site/year you want in the run.
#  (A single-site run is just a one-row table — no other code changes needed.)
wanted <- tibble::tribble(
  ~site, ~year,
  "CB", 19L,
  "TB", 19L,
  "TM", 19L,
  "GT", 19L,
  "HS", 19L,
  "N3", 19L
)

#"HS", 19L,
#"L1", 19L,
#"N3", 18L,
#"RS", 18L,
#"TB", 18L,
#"TM", 19L

env_long_test <- env_long %>%
  semi_join(wanted, by = c("site", "year")) %>%
  mutate(sitekey = make_key(site, year))

# the distinct site/year combos actually present after filtering (data-driven —
# guards against `wanted` containing a row that's missing from env_long)
site_years <- env_long_test %>%
  distinct(site, year) %>%
  mutate(sitekey = make_key(site, year)) %>%
  arrange(sitekey)

message("Run covers ", nrow(site_years), " site/year combos: ",
        paste(site_years$sitekey, collapse = ", "))

# ══════════════════════════════════════════════════════════════════════════════
#  RUN SUMMARY — print all active toggles before the slow work begins
#  Glance at this BEFORE the dem_store build kicks off (each prepare_dem is
#  expensive; better to catch a wrong toggle now than after 6 ambient_shades).
# ══════════════════════════════════════════════════════════════════════════════
cat("\n",
    "══════════════════════════════════════════════════════════════════\n",
    " RUN SUMMARY\n",
    "══════════════════════════════════════════════════════════════════\n",
    sep = "")
cat(sprintf("  Sites (%d)              : %s\n",
            nrow(site_years), paste(site_years$sitekey, collapse = ", ")))

cat("\n  ── Kd / attenuation ──\n")
cat(sprintf("  KD METHOD                : %s\n", kd_method_label()))
cat(sprintf("  KDPAR_COL                : %s\n", KDPAR_COL))
cat(sprintf("  USE_CONSTANT_KD          : %s\n", USE_CONSTANT_KD))
cat(sprintf("  USE_SPECTRAL_KD          : %s\n", USE_SPECTRAL_KD))
cat(sprintf("  USE_REFERENCE_KD         : %s\n", USE_REFERENCE_KD))
if (USE_REFERENCE_KD) {
  cat(sprintf("    N_REF_LOGGERS          : %d\n", N_REF_LOGGERS))
  for (sk in names(REFERENCE_LOGGER))
    cat(sprintf("    %-22s : %s\n", sk,
                paste(REFERENCE_LOGGER[[sk]], collapse = ", ")))
}
cat(sprintf("  KD_CONSTANT (fallback)   : %s\n", KD_CONSTANT))
cat("  Per-site KdPAR in use    :\n")
in_run <- intersect(site_years$sitekey, names(KDPAR_BY_KEY))
for (sk in in_run)
  cat(sprintf("    %-6s : %s\n", sk, round(KDPAR_BY_KEY[[sk]], 3)))
missing_kd <- setdiff(site_years$sitekey, names(KDPAR_BY_KEY))
if (length(missing_kd) > 0)
  cat(sprintf("    (no satellite Kd for: %s — will use fallback %.2f)\n",
              paste(missing_kd, collapse = ", "), KD_CONSTANT))

cat("\n  ── Direct/ambient weighting ──\n")
cat(sprintf("  USE_AMBIENT              : %s\n", USE_AMBIENT))
cat(sprintf("  USE_ALT_DEPENDENT_WEIGHTS: %s\n", USE_ALT_DEPENDENT_WEIGHTS))
cat(sprintf("  USE_DEPTH_DEP_WEIGHTS    : %s\n", USE_DEPTH_DEPENDENT_WEIGHTS))
cat(sprintf("  W_AMBIENT_MIN            : %s\n", W_AMBIENT_MIN))
cat(sprintf("  W_AMBIENT_MAX            : %s\n", W_AMBIENT_MAX))
cat(sprintf("  W_AMBIENT_P              : %s\n", W_AMBIENT_P))
cat(sprintf("  TAU_SCALE                : %s\n", TAU_SCALE))
cat(sprintf("  ALT_EXP                  : %s (>1 amplifies, <1 damps)\n", ALT_EXP))
cat(sprintf("  DEPTH_EXP                : %s (>1 amplifies, <1 damps)\n", DEPTH_EXP))
cat(sprintf("  WEIGHTS                  : %s\n", weight_label()))

cat("\n  ── DEM / raytrace geometry ──\n")
cat(sprintf("  REEF_RADIUS_M            : %s m\n", REEF_RADIUS_M))
cat(sprintf("  REEF_RADIUS_PX_TARGET    : %s px\n", REEF_RADIUS_PX_TARGET))
cat(sprintf("  MIN_DOWNSAMPLE           : %s\n", MIN_DOWNSAMPLE))
cat(sprintf("  MAXSEARCH_M              : %s m (direct beam)\n", MAXSEARCH_M))
cat(sprintf("  AMBIENT_MAXSEARCH_M      : %s m\n", AMBIENT_MAXSEARCH_M))
cat(sprintf("  AMBIENT_SUNBREAKS        : %s\n", AMBIENT_SUNBREAKS))
cat(sprintf("  BUFFER_M (logger extract): %s m\n", BUFFER_M))


cat("\n  ── Time ──\n")
cat(sprintf("  TIMESTEP_INTERVAL        : %s\n", TIMESTEP_INTERVAL))

cat("══════════════════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════════════════════════
#  LOAD ALL SPATIAL DATA ONCE  →  dem_store
#  This block runs prepare_dem (which includes ambient_shade — the slow step)
#  exactly ONCE per site/year. Every subsequent plot pulls from dem_store
#  instead of recomputing.
# ══════════════════════════════════════════════════════════════════════════════
dem_store <- list()
for (i in seq_len(nrow(site_years))) {
  s   <- site_years$site[i]
  y   <- site_years$year[i]
  key <- site_years$sitekey[i]
  
  message("\n── Loading spatial data for ", key, " ──")
  d   <- prepare_dem(s, y)
  hs  <- load_hobo_locations(s, y, geom)
  hs  <- extract_hobo_dem(d$rast, hs)
  
  dem_store[[key]] <- list(
    site       = s,
    year       = y,
    dem        = d,
    hobo_sf    = hs,
    hobo_utm   = st_transform(hs, UTM_CRS),
    hobo_wgs84 = st_transform(hs, "EPSG:4326")
  )
}


# ── OPTIONAL per-site DEM verification (histogram + ambient map), looped ───────
#  Set to FALSE to skip once you trust the DEMs. Useful when adding new sites
#  to spot-check that the reef mask & ambient layer look right.
VERIFY_DEMS <- FALSE
if (VERIFY_DEMS) {
  for (key in names(dem_store)) {
    e  <- dem_store[[key]]
    dm <- e$dem
    
    # distinguish real reef pixels from the "fill" border pixels (which are
    # set to reef_min - 100 m inside prepare_dem; well below any real value)
    fill_val    <- min(dm$matrix, na.rm = TRUE)
    real_pixels <- dm$matrix[dm$matrix != fill_val]
    
    dev.new()
    hist(real_pixels, breaks = 100,
         main = paste0(key, " DEM elevation (real pixels only)"),
         xlab = "Elevation (m)", col = "grey80", border = "white")
    
    cat("\n", key, " elevation quantiles:\n", sep = "")
    print(quantile(real_pixels, probs = c(0.01, 0.25, 0.50, 0.75, 0.99), na.rm = TRUE))
    cat("Total relief:", round(diff(range(real_pixels, na.rm = TRUE)), 3), "m\n")
    
    amb_rast_check <- rast(t(dm$ambient))
    ext(amb_rast_check) <- ext(dm$rast)
    crs(amb_rast_check) <- crs(dm$rast)
    amb_rast_check[!dm$reef_mask] <- NA
    dev.new()
    plot(amb_rast_check,
         main = paste0(key, " ambient occlusion (1=open sky, 0=occluded)"),
         col  = hcl.colors(100, "viridis"))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
#  RUN PIPELINE  (loops site/years internally; returns one pooled table)
# ══════════════════════════════════════════════════════════════════════════════
output <- run_pipeline(env_long_test, geom, dem_store) %>%
  mutate(sitekey = make_key(site, year))

# ══════════════════════════════════════════════════════════════════════════════
#  VALIDATION JOIN
#  Keyed on datetime + logger + SITEKEY so loggers can't collide across sites.
#  After the join, validation has one row per (sitekey, logger, hour) with
#  BOTH the prediction (pred_par, eff_light, shadow_fraction, ambient_fraction)
#  AND the measurement (light_mol). The filter drops hours where either is
#  missing (e.g. raytrace skipped at low sun, or sensor offline).
# ══════════════════════════════════════════════════════════════════════════════
validation <- output %>%
  left_join(
    env_long_test %>%
      mutate(hoboID_env = as.character(hoboID),
             datetime_pos = floor_date(datetime_pos, TIMESTEP_INTERVAL)) %>%
      group_by(datetime_pos, hoboID_env, sitekey) %>%
      summarise(light_mol = mean(light_mol, na.rm = TRUE), .groups = "drop"),
    by = c("datetime_pos", "hoboID_env", "sitekey")
  ) %>%
  filter(!is.na(pred_par), !is.na(light_mol))

# Drop reference loggers from validation when reference-Kd is on (they're
# INPUTS to Kd(t), so validating against them would be circular). The expansion
# step here is what required REFERENCE_LOGGER to be a list-of-vectors: a flat
# named vector would only drop one logger per site.
if (USE_REFERENCE_KD) {
  ref_drop <- tibble::tibble(
    sitekey    = rep(names(REFERENCE_LOGGER), lengths(REFERENCE_LOGGER)),
    hoboID_env = unlist(REFERENCE_LOGGER, use.names = FALSE)
  )
  validation <- validation %>% anti_join(ref_drop, by = c("sitekey", "hoboID_env"))
  message("Validation rows after dropping ", nrow(ref_drop),
          " reference logger(s): ", nrow(validation))
}

message("Validation rows: ", nrow(validation),
        " across ", length(unique(validation$sitekey)), " site/years")


# ══════════════════════════════════════════════════════════════════════════════
#  REFERENCE-LOGGER FINDER  (geometric openness + midday agreement)
#  ──────────────────────────────────────────────────────────────────────────
#  Picks N representative loggers per site to invert Beer-Lambert for Kd(t).
#  Score = weighted blend of TWO families of evidence:
#####
#  GEOMETRIC OPENNESS (the model says it's open) — from shadow & ambient
#    • shadow_mean : average openness across all daytime hours
#    • shadow_min  : worst-case openness (penalises loggers shaded part of day)
#    • shadow_sd   : steadiness (low variance = stable signal); we want loggers
#                    that are open ALL day, not just at noon and shaded at the
#                    ends — those give a biased Kd
#    • ambient     : open sky view (less ambient occlusion = sees full sky)
#
#  MIDDAY AGREEMENT (the data agrees the model is right about openness)
#    • pred/meas ratio at midday hours; closer to 1.0 = better
#    • A logger the MODEL thinks is open but the DATA says is dim is exactly
#      the WRONG reference — inverting Beer-Lambert there charges the missing
#      light to attenuation and inflates Kd. Adding this axis penalises that.
#    • Distance is in log space, so ratio 2.0 and ratio 0.5 score equally badly
#      (symmetric in over- vs under-prediction).
#
#  GATES (hard cutoffs — score forced to 0 if any fail):
#    • depth_ok   : finite positive mean depth with stable hour-to-hour values
#    • hours_ok   : enough daytime coverage to trust the summary
#    • midday_ok  : at least 2 midday hours with both pred & meas finite
#
#  WHY THIS REQUIRES `validation`, NOT `output`
#  ─────────────────────────────────────────────
#  We need light_mol to compute pred/meas. That's only in validation. The cost
#  is that the recommendation now depends on the run's Kd method (since pred_par
#  does). MITIGATION: run ONCE with USE_REFERENCE_KD = FALSE to anchor the
#  recommendation against satellite-Kd predictions, paste the recommendation
#  into REFERENCE_LOGGER, flip back to TRUE, and leave the IDs fixed.
#
#  NB: still a RECOMMENDATION — verify the top picks against the shadow maps
#  (is it actually on open sand, not perched on a coral head?).
# ══════════════════════════════════════════════════════════════════════════════
find_reference_loggers <- function(val,
                                   n_top         = N_REF_LOGGERS,
                                   depth_sd_max  = 3,    # max hour-to-hour depth sd (m)
                                   min_hours     = 5,    # min daytime hours to qualify
                                   midday_hours  = 10:14,# local hours used for agreement
                                   # geometric-openness axis weights
                                   w_shadow_mean = 0.4, # weight: average openness to sun
                                   w_shadow_min  = 0.15, # weight: worst-case openness
                                   w_shadow_stab = 0.10, # weight: temporal steadiness
                                   w_ambient     = 0.10, # weight: openness to sky
                                   # data-agreement axis weight
                                   w_agreement   = 0.25) {
  stopifnot(abs(w_shadow_mean + w_shadow_min + w_shadow_stab +
                  w_ambient + w_agreement - 1) < 1e-9)
  
  # per-logger summary: openness statistics + midday pred/meas ratio
  scored <- val %>%
    filter(!is.na(shadow_fraction)) %>%      # daytime, raytraced rows only
    mutate(hour = as.integer(format(datetime_pos, "%H"))) %>%
    group_by(sitekey, hoboID_env) %>%
    summarise(
      shadow_mean  = mean(shadow_fraction,  na.rm = TRUE),
      shadow_min   = min(shadow_fraction,   na.rm = TRUE),  # worst hour of the day
      shadow_sd    = sd(shadow_fraction,    na.rm = TRUE),  # day-cycle variability
      ambient_mean = mean(ambient_fraction, na.rm = TRUE),
      depth_mean   = mean(water_depth_m,    na.rm = TRUE),
      depth_sd     = sd(water_depth_m,      na.rm = TRUE),
      # midday-only pred/meas ratio. Sum (not ratio of means) is a proxy for the
      # daily light dose at midday; identical to ratio-of-means but composes
      # better when some hours are NA on one side.
      midday_pred  = sum(pred_par[hour %in% midday_hours],  na.rm = TRUE),
      midday_meas  = sum(light_mol[hour %in% midday_hours], na.rm = TRUE),
      n_midday     = sum(hour %in% midday_hours & !is.na(pred_par) &
                           !is.na(light_mol)),
      n_hours      = dplyr::n_distinct(datetime_pos),
      .groups = "drop"
    ) %>%
    mutate(
      midday_ratio = ifelse(midday_meas > 0 & midday_pred > 0,
                            midday_pred / midday_meas, NA_real_),
      # log-distance from 1.0 — symmetric (ratio 2 and ratio 0.5 score equally).
      # negative = worse, max = 0 (perfect agreement). norm01 below will flip
      # this so higher = better, matching the other axes.
      agreement_raw = ifelse(is.finite(midday_ratio) & midday_ratio > 0,
                             -abs(log(midday_ratio)), NA_real_),
      # HARD GATES
      depth_ok  = is.finite(depth_mean) & depth_mean > 0 &
        is.finite(depth_sd)   & depth_sd < depth_sd_max,
      hours_ok  = n_hours >= min_hours,
      midday_ok = n_midday >= 2 & is.finite(agreement_raw),
      gate_ok   = depth_ok & hours_ok & midday_ok
    )
  
  # Normalise each axis to 0–1 WITHIN site so weights mean the same thing
  # everywhere (each site's own best logger sits at ~1 per axis; if Site A is
  # less open than Site B overall, we still pick Site A's relative best).
  # shadow_sd is inverted: low variance = steady = good -> high stability_norm.
  norm01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    o <- (x - r[1]) / (diff(r) + 1e-9)
    o[!is.finite(o)] <- 0   # non-finite (e.g. sd of 1 obs) -> worst, never NA
    o
  }
  scored <- scored %>%
    group_by(sitekey) %>%
    mutate(
      shadow_mean_norm = norm01(shadow_mean),
      shadow_min_norm  = norm01(shadow_min),
      stability_norm   = 1 - norm01(shadow_sd),
      ambient_norm     = norm01(ambient_mean),
      agreement_norm   = norm01(agreement_raw),   # best site logger -> 1, worst -> 0
      open_score = ifelse(
        gate_ok,
        w_shadow_mean * shadow_mean_norm +
          w_shadow_min  * shadow_min_norm  +
          w_shadow_stab * stability_norm   +
          w_ambient     * ambient_norm     +
          w_agreement   * agreement_norm,
        0)
    ) %>%
    arrange(sitekey, desc(open_score)) %>%
    ungroup()
  
  # ranked table per site (so you can eyeball the alternatives and not just the winner)
  for (k in sort(unique(scored$sitekey))) {
    cat("\n── Reference candidates for ", k, " (ranked) ──\n", sep = "")
    sub <- scored %>% filter(sitekey == k) %>%
      transmute(hoboID_env,
                shadow       = round(shadow_mean, 3),
                shadow_min   = round(shadow_min,  3),
                ambient      = round(ambient_mean, 3),
                depth_m      = round(depth_mean,  2),
                midday_ratio = round(midday_ratio, 2),
                n_midday,
                n_hours,
                gate_ok,
                score        = round(open_score, 3))
    print(as.data.frame(sub), row.names = FALSE)
  }
  
  # top N per site (must pass gates and have positive score)
  best <- scored %>%
    filter(open_score > 0) %>%
    group_by(sitekey) %>%
    slice_max(open_score, n = n_top, with_ties = FALSE) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    select(sitekey, rank, hoboID_env, shadow_mean, midday_ratio,
           depth_mean, n_hours, open_score)
  
  # flag any site that produced NO qualifying logger (so it's not silently missing)
  missing_sites <- setdiff(sort(unique(scored$sitekey)), best$sitekey)
  if (length(missing_sites) > 0)
    message("  [ref-finder] NO qualifying reference logger for: ",
            paste(missing_sites, collapse = ", "),
            " (relax depth_sd_max / min_hours, or accept satellite Kd there)")
  
  cat("\n── TOP ", n_top,
      " reference loggers per site (openness + midday agreement) ──\n", sep = "")
  print(as.data.frame(best %>%
                        transmute(sitekey, rank, hoboID_env,
                                  shadow       = round(shadow_mean, 3),
                                  midday_ratio = round(midday_ratio, 2),
                                  depth_m      = round(depth_mean, 2),
                                  n_hours,
                                  score        = round(open_score, 3))),
        row.names = FALSE)
  
  # ready-to-paste list for REFERENCE_LOGGER
  cat("\n# Paste into REFERENCE_LOGGER (top ", n_top, " per site):\n", sep = "")
  cat("REFERENCE_LOGGER <- list(\n")
  site_keys <- unique(best$sitekey)
  for (i in seq_along(site_keys)) {
    k   <- site_keys[i]
    ids <- best$hoboID_env[best$sitekey == k]
    cat(sprintf('  "%s" = c(%s)%s\n',
                k,
                paste0('"', ids, '"', collapse = ", "),
                if (i < length(site_keys)) "," else ""))
  }
  cat(")\n")
  
  invisible(best)
}

# Score from validation (needs light_mol). See bootstrap workflow in the
# REFERENCE_LOGGER block above.
ref_candidates <- find_reference_loggers(validation)#start
#######


# ── RUN TAG (set once; used in every saved PNG filename + title) ──────────────
RUN_TAG     <- format(Sys.time(), "%Y%m%d_%H%M%S")
RUNS_ROOT   <- "runs"                              # parent folder; holds metrics_log.csv
RUNS_DIR    <- file.path(RUNS_ROOT, RUN_TAG)       # this run's own folder
dir.create(RUNS_DIR, showWarnings = FALSE, recursive = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
#  PER-SITE FIGURE 1 — measured vs predicted PAR per logger
#  One figure PER site/year, saved to file. Loops over dem_store keys.
#  Useful to spot specific loggers where the prediction badly diverges from the
#  measurement — they usually have one of: wrong placement (lon/lat off),
#  wrong dpHB, or unmodelled overhang (3-D cavity the top-down DEM flattens).
# ══════════════════════════════════════════════════════════════════════════════
plot_par_panels_one_site <- function(key, val_all) {
  val <- val_all %>% filter(sitekey == key)
  if (nrow(val) == 0) { message("  (no validation rows for ", key, ")"); return(invisible()) }
  
  date_range  <- format(range(val$datetime_pos), "%Y-%m-%d %H:%M")
  out_px      <- max(dim(dem_store[[key]]$dem$matrix))
  # Kd label: per-hour if reference-Kd, satellite value if not, fallback otherwise
  kd_disp <- if (USE_REFERENCE_KD) "per-hour"
  else ifelse(key %in% names(KDPAR_BY_KEY), round(KDPAR_BY_KEY[[key]],3), KD_CONSTANT)
  plot_title  <- paste0(
    key, " | ", RUN_TAG, " | ", kd_method_label(), " | ",
    date_range[1], " to ", date_range[2], "\n",
    "Kd=", kd_disp, " | zscale=", round(dem_store[[key]]$dem$zscale, 4), " (UTM)",
    " | interval=", TIMESTEP_INTERVAL, " | ~", out_px, "px matrix",
    " | reef=", REEF_RADIUS_M, "m (~", REEF_RADIUS_PX_TARGET, "px target)",
    ", maxsearch=", MAXSEARCH_M, "m | ", weight_label()
  )
  
  # shared y-axis range so panels are visually comparable
  ylim_shared <- range(c(val$pred_par, val$light_mol), na.rm = TRUE)
  hobos       <- sort(unique(val$hoboID_env))
  n_hobos     <- length(hobos)
  n_cols      <- 5
  n_rows      <- ceiling(n_hobos / n_cols)
  legend_cell <- n_hobos + 2
  plot_cells  <- c(2:(n_hobos + 1), rep(0, n_cols * n_rows - n_hobos))
  
  layout_mat <- rbind(
    rep(1, n_cols),
    matrix(plot_cells, nrow = n_rows, ncol = n_cols, byrow = TRUE),
    rep(legend_cell, n_cols)
  )
  
  png(file.path(RUNS_DIR, paste0("par_panels_", key, "_", RUN_TAG, ".png")),
      width = n_cols * 360, height = n_rows * 300 + 220, res = 130)
  layout(layout_mat, heights = c(0.8, rep(3, n_rows), 0.8))
  
  par(mar = c(0, 0, 0, 0)); plot.new()
  text(0.5, 0.5, plot_title, cex = 1.1, font = 2, adj = 0.5)
  
  for (j in seq_along(hobos)) {
    par(mar = c(2, 2, 1.5, 0.5))
    id  <- hobos[j]
    sub <- val[val$hoboID_env == id, ]
    sub <- sub[order(sub$datetime_pos), ]
    is_left <- (j - 1) %% n_cols == 0
    
    plot(sub$datetime_pos, sub$light_mol,
         type = "l", col = "blue", lwd = 1.5, ylim = ylim_shared,
         xlab = "", ylab = "", main = id, cex.main = 0.85,
         xaxt = "n", yaxt = if (is_left) "s" else "n")
    axis.POSIXct(1, sub$datetime_pos, format = "%H", cex.axis = 0.6)
    lines(sub$datetime_pos, sub$pred_par, col = "red", lwd = 1.5, lty = 2)
  }
  
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = c("Measured", "Predicted"),
         col = c("blue", "red"), lty = c(1, 2), lwd = 2,
         cex = 1.0, bty = "n", horiz = TRUE)
  layout(1); dev.off()
  message("  saved par_panels_", key, ".png")
}

for (key in names(dem_store)) plot_par_panels_one_site(key, validation)


# ══════════════════════════════════════════════════════════════════════════════
#  PER-SITE FIGURE 2 — shadow maps + per-timestep bar charts
#  One figure PER site/year, saved to file. Reuses cached DEM/hobo from dem_store.
#  Useful to confirm a logger sits on the shadow map where you expect AND to
#  see how the predictions compare to measurements at each specific hour.
# ══════════════════════════════════════════════════════════════════════════════
plot_shadow_multipanel_one_site <- function(key, output_all, val_all) {
  e          <- dem_store[[key]]
  dem        <- e$dem
  hobo_sf    <- e$hobo_sf
  hobo_wgs84 <- e$hobo_wgs84
  
  out <- output_all %>% filter(sitekey == key)
  val <- val_all    %>% filter(sitekey == key)
  if (nrow(out) == 0) { message("  (no output rows for ", key, ")"); return(invisible()) }
  
  # one timestep per row -> shadow map + bar chart of pred vs meas
  sun_steps <- out %>% distinct(datetime_pos, sun_altitude, sun_azimuth) %>%
    arrange(datetime_pos)
  n_steps <- nrow(sun_steps)
  n_cols  <- 4
  n_rows  <- ceiling(n_steps / n_cols)
  
  # per-DEM maxsearch in pixels — same conversion as make_shadow, so the
  # plotted shadow maps use IDENTICAL ray distances to the pipeline output.
  maxsearch_px_plot <- max(1, round(MAXSEARCH_M / dem$zscale))
  
  draw <- function() {
    layout_mat <- matrix(0, nrow = n_rows * 2 + 2, ncol = n_cols)
    layout_mat[1, ] <- 1
    panel <- 2
    for (r in seq_len(n_rows)) for (c in seq_len(n_cols)) {
      step_i <- (r - 1) * n_cols + c
      if (step_i <= n_steps) {
        layout_mat[r * 2,     c] <- panel
        layout_mat[r * 2 + 1, c] <- panel + 1
        panel <- panel + 2
      }
    }
    legend_cell <- panel
    layout_mat[n_rows * 2 + 2, ] <- legend_cell
    layout(layout_mat, heights = c(0.6, rep(c(5, 1.5), n_rows), 1.5))
    
    par(mar = c(0, 0, 0, 0)); plot.new()
    text(0.5, 0.5,
         paste0(key, "  |  ", RUN_TAG, "  |  ", kd_method_label(),
                "  |  Kd=", if (USE_REFERENCE_KD) "per-hour"
                else ifelse(key %in% names(KDPAR_BY_KEY), round(KDPAR_BY_KEY[[key]],3), KD_CONSTANT),
                "  | zscale=", round(dem$zscale, 4), " (UTM)",
                "  |  ~", max(dim(dem$matrix)), "px matrix",
                "  |  reef=", REEF_RADIUS_M, "m (~",
                round(REEF_RADIUS_M / dem$zscale), "px) ",
                "maxsearch=", MAXSEARCH_M, "m (", maxsearch_px_plot, "px)",
                "  |  ", weight_label()),
         cex = 1.2, font = 2, adj = 0.5)
    
    for (i in seq_len(n_steps)) {
      ts      <- sun_steps[i, ]
      alt_deg <- ts$sun_altitude * 180 / pi
      az_deg  <- ts$sun_azimuth  * 180 / pi
      
      par(mar = c(0, 1, 5, 1))
      # re-raytrace for the plot, using the SAME physical maxsearch as the
      # pipeline (converted to pixels above via MAXSEARCH_M / zscale)
      shadow <- ray_shade(dem$matrix, sunaltitude = alt_deg, sunangle = az_deg,
                          zscale = dem$zscale, maxsearch = maxsearch_px_plot)
      shadow[!dem$reef_mask] <- NA
      shadow_rast <- rast(shadow)
      ext(shadow_rast) <- ext(dem$rast)
      crs(shadow_rast) <- crs(dem$rast)
      # trim NA padding so the disc fills the panel even when the raw DEM
      # was much wider than the 6.5 m reef radius (e.g. CB19 raw extent ~23 m
      # vs CB18 ~13 m). Physics is identical between sites; this is purely
      # cosmetic so plots are visually comparable.
      shadow_rast <- terra::trim(shadow_rast)
      
      plot(project(shadow_rast, "EPSG:4326"),
           main = paste0(format(ts$datetime_pos, "%Y-%m-%d %H:%M"),
                         "\nalt=", round(alt_deg, 1), "° az=", round(az_deg, 1), "°"),
           legend = FALSE, axes = FALSE, cex.main = 1.0)
      plot(vect(hobo_wgs84), add = TRUE, col = "red", bg = NA,
           cex = 1.5, pch = 21, lwd = 2)
      text(st_coordinates(hobo_wgs84)[, 1], st_coordinates(hobo_wgs84)[, 2],
           labels = hobo_sf$hoboID, col = "yellow", cex = 0.9, font = 2, pos = 3)
      
      par(mar = c(4, 3, 0.5, 1))
      step_data <- val %>% filter(datetime_pos == ts$datetime_pos) %>%
        arrange(hoboID_env)
      if (nrow(step_data) > 0) {
        bar_mat <- rbind(step_data$pred_par, step_data$light_mol)
        bar_mat[is.na(bar_mat)] <- 0
        barplot(bar_mat, beside = TRUE, col = c("red", "blue"),
                names.arg = step_data$hoboID_env, cex.names = 0.8,
                cex.axis = 0.9, las = 2,
                ylim = c(0, max(c(val$pred_par, val$light_mol), na.rm = TRUE)),
                border = NA)
      } else plot.new()
    }
    
    # legend / colour key strip at the bottom
    par(mar = c(2, 4, 2, 4)); plot.new()
    n_grad <- 100; grad <- hcl.colors(n_grad, "viridis")
    x_left <- 0.05; x_right <- 0.55; y_bot <- 0.3; y_top <- 0.65
    x_steps <- seq(x_left, x_right, length.out = n_grad + 1)
    for (j in seq_len(n_grad))
      rect(x_steps[j], y_bot, x_steps[j + 1], y_top, col = grad[j], border = NA)
    rect(x_left, y_bot, x_right, y_top, border = "black", lwd = 0.5)
    text(x_left,  y_bot - 0.15, "0\n(full shadow)", cex = 1.1, adj = 0.5)
    text(x_right, y_bot - 0.15, "1\n(full sun)",    cex = 1.1, adj = 0.5)
    text((x_left + x_right) / 2, y_top + 0.15, "Shadow fraction", cex = 1.3, font = 2)
    legend(0.65, 0.8, legend = c("Predicted PAR", "Measured PAR"),
           fill = c("red", "blue"), border = NA, bty = "n", cex = 1.2,
           title = "Bar chart key")
  }
  
png(file.path(RUNS_DIR, paste0("shadow_multipanel_", key, "_", RUN_TAG, ".png")),
      width = n_cols * 800, height = n_rows * 1100 + 200, res = 150)
  draw(); dev.off()
  message("  saved shadow_multipanel_", key, ".png")
}

for (key in names(dem_store)) plot_shadow_multipanel_one_site(key, output, validation)



# ══════════════════════════════════════════════════════════════════════════════
#  POOLED-PLOT SITE FILTER  (applies to BOTH pooled figures below)
#  ──────────────────────────────────────────────────────────────────────────
#  Default "all" pools every site/year in the run. To show only some, set this
#  to a vector of sitekeys, e.g.  POOL_SITES <- c("TB18")  or  c("CB18","CB19").
#  Keys must match the sitekey format ("TB18", "CB19", ...). The per-site PAR
#  panels and shadow maps above are unaffected — they always loop all sites.
# ══════════════════════════════════════════════════════════════════════════════
POOL_SITES <- "all"        # "all"  OR  e.g. c("TB18")  OR  c("CB18","CB19")

# resolve the filter to a concrete set of keys, with a guard for typos
pool_keys_all <- sort(unique(validation$sitekey))
pool_keys <- if (identical(POOL_SITES, "all")) pool_keys_all else POOL_SITES

bad_keys <- setdiff(pool_keys, pool_keys_all)
if (length(bad_keys) > 0)
  warning("POOL_SITES has keys not present in this run: ",
          paste(bad_keys, collapse = ", "),
          " | available: ", paste(pool_keys_all, collapse = ", "))
pool_keys <- intersect(pool_keys, pool_keys_all)
if (length(pool_keys) == 0)
  stop("POOL_SITES matched no site/years in the run. Available: ",
       paste(pool_keys_all, collapse = ", "))

message("Pooled plots include: ", paste(pool_keys, collapse = ", "))


# ══════════════════════════════════════════════════════════════════════════════
#  POOLED FIGURE 1 — DLI predicted vs measured
#  SHAPE = site/year, COLOUR = depth (or shadow), every sensor labelled.
#  This is the overall "how well does the model do on cumulative light" check.
# ══════════════════════════════════════════════════════════════════════════════
# ── COLOUR TOGGLE: "depth" or "shadow" ────────────────────────────────────────
COLOR_BY <- "depth"   # "depth" = water column (m); "shadow" = mean shadow fraction

# ── DLI UNIT / INTEGRATION-WINDOW FIX ─────────────────────────────────────────
#  OLD BUG (fixed): measured DLI summed the raw 5-min light_mol from
#  env_long_test (× 300 s) while predicted DLI summed hourly pred_par (× 3600 s).
#  The two therefore integrated over DIFFERENT row sets and DIFFERENT time bases:
#    • measured included every 5-min slot inside a qualifying HOUR, even slots
#      where pred_par was NA (raytrace skipped, depth<=0) — so a sensor could
#      accumulate measured light in hours it had no prediction for.
#    • the 300 vs 3600 factors were each internally correct for their own path,
#      but the MISMATCHED windows made measured-vs-predicted DLI not comparable.
#
#  FIX: integrate BOTH from `validation`, which already holds ONE hourly-mean
#  row per logger/hour for both light_mol and pred_par (the validation join
#  averaged light_mol over the same hour pred_par is built on). Both are mean
#  PAR over the hour in umol m^-2 s^-1, so BOTH multiply by 3600 s and divide by
#  1e6 -> mol m^-2. Identical rows in, identical window out: directly comparable.
#  Only hours present in `validation` (sun_alt>20 AND a finite prediction) count
#  on both sides, which is exactly what we want.
HOUR_SECONDS <- 3600        # validation rows are hourly means

dli <- validation %>%
  filter(sitekey %in% pool_keys) %>%
  group_by(sitekey, hoboID_env) %>%
  summarise(
    dli_measured  = sum(light_mol, na.rm = TRUE) * HOUR_SECONDS / 1e6,
    dli_predicted = sum(pred_par,  na.rm = TRUE) * HOUR_SECONDS / 1e6,
    depth_m       = mean(water_depth_m,    na.rm = TRUE),
    shadow_mean   = mean(shadow_fraction,  na.rm = TRUE),
    n_hours       = n(),
    .groups = "drop"
  ) %>%
  as.data.frame()

# SHAPE per site/year (pch_map maps sitekey -> a unique pch value)
keys    <- sort(unique(dli$sitekey))
pch_map <- setNames(seq_along(keys) %% 25, keys)

# ── COLOUR driven by the toggle ───────────────────────────────────────────────
if (COLOR_BY == "shadow") {
  col_vals <- dli$shadow_mean; col_lab <- "Shadow"; col_dig <- 2
} else {
  col_vals <- dli$depth_m;     col_lab <- "Depth (m)"; col_dig <- 1
}
col_rng <- range(col_vals, na.rm = TRUE)
col_pal <- hcl.colors(100, "viridis", rev = TRUE)
pt_col  <- col_pal[
  pmax(1, pmin(100, round((col_vals - col_rng[1]) / diff(col_rng) * 99) + 1))
]


png(file.path(RUNS_DIR, sprintf("pooled_dli_%s.png", RUN_TAG)),
    width = 8 * 130, height = 7 * 130, res = 130)
par(mar = c(5, 5, 4, 7), xpd = FALSE)
lim <- range(c(dli$dli_measured, dli$dli_predicted), na.rm = TRUE)

plot(dli$dli_measured, dli$dli_predicted,
     pch = pch_map[dli$sitekey], col = pt_col, cex = 1.5, lwd = 2,
     xlim = lim, ylim = lim,
     xlab = "Measured total light (mol m-2, deployment sum)",
     ylab = "Predicted total light (mol m-2, deployment sum)",
     main = paste0("Total summed daytime light measured vs predicted — ", length(keys),
                   " site/years (colour = ", col_lab, ") | ", kd_method_label(),
                   "\nsatellite KdPAR per-site: ",
                   paste(sprintf("%s=%.2f", names(KDPAR_BY_KEY)[match(keys, names(KDPAR_BY_KEY))],
                                 KDPAR_BY_KEY[match(keys, names(KDPAR_BY_KEY))]), collapse=", "),
                   if (USE_REFERENCE_KD) " (shown for ref; actual Kd is per-hour)" else "",
                   " | reef=", REEF_RADIUS_M, "m, maxsearch=", MAXSEARCH_M, "m",
                   " | ", weight_label()),
     cex.main = 0.8)
abline(0, 1, lty = 2, col = "red", lwd = 2)
text(dli$dli_measured, dli$dli_predicted, labels = dli$hoboID_env,
     cex = 0.6, pos = 3)

# shape legend (site/year)
par(xpd = TRUE)
legend(lim[2] + diff(lim) * 0.04, lim[2],
       legend = keys, pch = pch_map[keys], pt.cex = 1.4, pt.lwd = 2,
       bty = "n", title = "Site/year")

# colour strip (label + range follow the toggle)
strip_x <- lim[2] + diff(lim) * 0.04
strip_y <- seq(lim[1], lim[1] + diff(lim) * 0.45, length.out = 100)
for (k in 1:99)
  rect(strip_x, strip_y[k], strip_x + diff(lim) * 0.03, strip_y[k + 1],
       col = col_pal[k], border = NA)
text(strip_x, max(strip_y) + diff(lim) * 0.03, col_lab, cex = 0.8, font = 2, adj = 0)
text(strip_x + diff(lim) * 0.04, min(strip_y), round(col_rng[1], col_dig), cex = 0.7, adj = 0)
text(strip_x + diff(lim) * 0.04, max(strip_y), round(col_rng[2], col_dig), cex = 0.7, adj = 0)
par(xpd = FALSE)
dev.off()



# ══════════════════════════════════════════════════════════════════════════════
#  POOLED FIGURE 2 — observed vs predicted PAR at selected hours
#  Shows per-hour bias (does the model under-predict at 17h? over-predict at 12h?).
#  Dotted lines are slope-through-origin fits per hour, so a hour's points all
#  lying above its dotted line means the model over-predicts at that hour.
# ══════════════════════════════════════════════════════════════════════════════
time_df <- validation %>%
  filter(sitekey %in% pool_keys) %>%
  mutate(hour = as.integer(format(datetime_pos, "%H"))) %>%
  filter(hour %in% c(8, 10, 12, 15, 17))

if (nrow(time_df) > 0) {
  lim2 <- range(c(time_df$light_mol, time_df$pred_par), na.rm = TRUE)
  
  keys2    <- sort(unique(time_df$sitekey))
  pch_map2 <- setNames(seq_along(keys2) %% 25, keys2)
  
  # COLOUR = hour of day (the variable we're investigating)
  hour_levels <- c(8, 10, 12, 15, 17)
  hour_cols   <- c( "8"= "green",
                    "10" = "coral",
                    "12" = "#E69F00",   # midday  — amber
                    "15" = "#0072B2",   # mid-aft  — blue
                    "17" = "#CC79A7")   # late aft — magenta
  col2 <- hour_cols[as.character(time_df$hour)]
  
  
  png(file.path(RUNS_DIR, sprintf("pooled_perhour_%s.png", RUN_TAG)),
      width = 8 * 130, height = 7 * 130, res = 130)
  par(mar = c(5, 5, 4, 7))
  plot(time_df$light_mol, time_df$pred_par,
       pch = pch_map2[time_df$sitekey], col = col2, cex = 1.4, lwd = 2,
       xlim = lim2, ylim = lim2,
       xlab = "Observed PAR (µmol)", ylab = "Predicted PAR (µmol)",
       main = paste0("PAR at 8/10/12/15/17h — ", length(keys2),
                     " site/years (colour = hour) | ", kd_method_label(),
                     " | ", weight_label()),
       cex.main = 0.85)
  abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
  text(time_df$light_mol, time_df$pred_par, labels = time_df$hoboID_env,
       cex = 0.55, pos = 3)
  
  # per-hour bias guides: a fitted slope-through-origin for each hour, so you can
  # SEE whether a given hour systematically over- or under-predicts.
  for (h in hour_levels) {
    sub <- time_df[time_df$hour == h, ]
    if (nrow(sub) > 2) {
      b  <- coef(lm(pred_par ~ 0 + light_mol, data = sub))[1]
      xr <- range(sub$light_mol, na.rm = TRUE)        # only this hour's x-range
      segments(xr[1], b*xr[1], xr[2], b*xr[2],
               col = hour_cols[as.character(h)], lwd = 2, lty = 3)
    }
  }
  
  par(xpd = TRUE)
  # hour legend (colour) — the primary key
  legend(lim2[2] + diff(lim2) * 0.04, lim2[2],
         legend = paste0(hour_levels, ":00"),
         col = hour_cols[as.character(hour_levels)],
         pch = 19, pt.cex = 1.4, bty = "n", title = "Hour")
  # site legend (shape)
  legend(lim2[2] + diff(lim2) * 0.04, lim2[2] - diff(lim2) * 0.28,
         legend = keys2, pch = pch_map2[keys2], pt.cex = 1.3, pt.lwd = 2,
         bty = "n", title = "Site/year")
  par(xpd = FALSE)
  dev.off()
  
}

# ══════════════════════════════════════════════════════════════════════════════
#  RUN OUTPUTS — per-run manifest + metrics + accumulating sweep log
#  ──────────────────────────────────────────────────────────────────────────
#  Each run produces THREE artefacts in RUNS_DIR:
#    1. manifest_<RUN_TAG>.txt   — human-readable settings + metrics summary
#    2. metrics_<RUN_TAG>.csv    — full per-site + pooled metrics table
#    3. metrics_log.csv          — APPENDS one summary row per run (sweeps)
#
#  RUN_TAG is a timestamp so nothing overwrites between runs.
# ══════════════════════════════════════════════════════════════════════════════

# ── metrics helper ────────────────────────────────────────────────────────────
metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  n   <- length(obs)
  if (n < 2) return(data.frame(n=n, rmse=NA, mae=NA, bias=NA,
                               r2=NA, r2_1to1=NA, slope=NA, intercept=NA))
  res    <- pred - obs
  rmse   <- sqrt(mean(res^2))
  mae    <- mean(abs(res))
  bias   <- mean(res)                              # +ve = over-predicts
  ss_res <- sum(res^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2_1to1 <- 1 - ss_res / ss_tot                   # Nash-Sutcliffe vs 1:1
  fit    <- lm(pred ~ obs)
  r2     <- summary(fit)$r.squared                 # Pearson R² (best-fit line)
  data.frame(n=n, rmse=rmse, mae=mae, bias=bias,
             r2=r2, r2_1to1=r2_1to1,
             slope=coef(fit)[2], intercept=coef(fit)[1])
}

# ── (a) HOURLY metrics ────────────────────────────────────────────────────────
hourly_per_site <- validation %>%
  group_by(sitekey) %>%
  group_modify(~ metrics(.x$light_mol, .x$pred_par)) %>%
  ungroup() %>% mutate(level = "hourly_per_site")

hourly_overall <- metrics(validation$light_mol, validation$pred_par) %>%
  mutate(sitekey = "ALL", level = "hourly_overall")

# ── (b) DLI metrics (deployment sums, mol m^-2) ───────────────────────────────
dli_per_site <- dli %>%
  group_by(sitekey) %>%
  group_modify(~ metrics(.x$dli_measured, .x$dli_predicted)) %>%
  ungroup() %>% mutate(level = "dli_per_site")

dli_overall <- metrics(dli$dli_measured, dli$dli_predicted) %>%
  mutate(sitekey = "ALL", level = "dli_overall")

metrics_tbl <- bind_rows(hourly_overall, hourly_per_site,
                         dli_overall,    dli_per_site) %>%
  select(level, sitekey, n, rmse, mae, bias, r2, r2_1to1, slope, intercept)

# ── print to console ──────────────────────────────────────────────────────────
cat("\n================ VALIDATION METRICS ================\n")
cat(sprintf("Run tag: %s\n", RUN_TAG))
cat(sprintf("Kd method: %s | weights: %s\n", kd_method_label(), weight_label()))
cat(sprintf("REEF_RADIUS_M=%.1f  MAXSEARCH_M=%.1f  AMBIENT_MAXSEARCH_M=%.2f\n",
            REEF_RADIUS_M, MAXSEARCH_M, AMBIENT_MAXSEARCH_M))
print(metrics_tbl, row.names = FALSE, digits = 3)

# ── save per-run metrics table ────────────────────────────────────────────────
write.csv(metrics_tbl,
          file.path(RUNS_DIR, sprintf("metrics_%s.csv", RUN_TAG)),
          row.names = FALSE)

# ── HUMAN-READABLE PER-RUN MANIFEST ───────────────────────────────────────────
manifest_path <- file.path(RUNS_DIR, sprintf("manifest_%s.txt", RUN_TAG))
sink(manifest_path)

cat("RUN TAG : ", RUN_TAG, "\n", sep = "")
cat("RUN DATE: ", format(Sys.time()), "\n\n", sep = "")

cat("SITES:\n")
cat("  ", paste(sort(unique(validation$sitekey)), collapse = ", "), "\n\n", sep = "")

cat("KD METHOD : ", kd_method_label(), "\n", sep = "")
cat("WEIGHTING : ", weight_label(),    "\n\n", sep = "")

cat("SETTINGS:\n")
cat("  USE_CONSTANT_KD           =", USE_CONSTANT_KD, "\n")
cat("  USE_REFERENCE_KD          =", USE_REFERENCE_KD, "\n")
cat("  USE_SPECTRAL_KD           =", USE_SPECTRAL_KD, "\n")
cat("  USE_AMBIENT               =", USE_AMBIENT, "\n")
cat("  USE_ALT_DEPENDENT_WEIGHTS =", USE_ALT_DEPENDENT_WEIGHTS, "\n\n")
cat("  USE_DEPTH_DEPENDENT_WEIGHTS =", USE_DEPTH_DEPENDENT_WEIGHTS, "\n\n")
cat("  TAU_SCALE                 =", TAU_SCALE, "\n\n")
cat("  ALT_EXP                   =", ALT_EXP, "\n")
cat("  DEPTH_EXP                 =", DEPTH_EXP, "\n")

cat("  REEF_RADIUS_M             =", REEF_RADIUS_M, "\n")
cat("  REEF_RADIUS_PX_TARGET     =", REEF_RADIUS_PX_TARGET, "\n")
cat("  MIN_DOWNSAMPLE            =", MIN_DOWNSAMPLE, "\n\n")

cat("  MAXSEARCH_M               =", MAXSEARCH_M, "\n")
cat("  AMBIENT_MAXSEARCH_M       =", AMBIENT_MAXSEARCH_M, "\n")
cat("  AMBIENT_SUNBREAKS         =", paste(AMBIENT_SUNBREAKS, collapse = ", "), "\n\n")

cat("  KDPAR_COL                 =", KDPAR_COL, "\n")
cat("  KD_CONSTANT               =", KD_CONSTANT, "\n")

cat("  W_DIRECT                  =", W_DIRECT, "\n")
cat("  W_AMBIENT                 =", W_AMBIENT, "\n")
cat("  W_AMBIENT_MIN             =", W_AMBIENT_MIN, "\n")
cat("  W_AMBIENT_MAX             =", W_AMBIENT_MAX, "\n")
cat("  W_AMBIENT_P               =", W_AMBIENT_P, "\n\n")

cat("  TIMESTEP_INTERVAL         =", TIMESTEP_INTERVAL, "\n")
cat("  BUFFER_M                  =", BUFFER_M, "\n")
cat("  KDPAR_COL                 =", KDPAR_COL, "\n")

if (USE_REFERENCE_KD) {
  cat("\nREFERENCE LOGGERS:\n")
  for (site in names(REFERENCE_LOGGER)) {
    cat("  ", site, ": ",
        paste(REFERENCE_LOGGER[[site]], collapse = ", "),
        "\n", sep = "")
  }
}

cat("\n================ VALIDATION METRICS ================\n")
print(metrics_tbl, row.names = FALSE, digits = 3)

sink()
message("Saved ", manifest_path)

# ── APPEND one-row summary to the accumulating sweep log ──────────────────────
log_row <- data.frame(
  run_tag             = RUN_TAG,
  sites               = paste(sort(unique(validation$sitekey)), collapse = ";"),
  kd_method           = kd_method_label(),
  kdpar_col           = KDPAR_COL,
  use_spectral_kd     = USE_SPECTRAL_KD,
  use_reference_kd    = USE_REFERENCE_KD,
  use_ambient         = USE_AMBIENT,
  use_alt_dep_weights = USE_ALT_DEPENDENT_WEIGHTS,
  use_depth_dep_weights = USE_DEPTH_DEPENDENT_WEIGHTS,
  tau_scale             = TAU_SCALE,
  alt_exp               = ALT_EXP,
  depth_exp            = DEPTH_EXP,
  w_direct            = W_DIRECT,
  w_ambient           = W_AMBIENT,
  w_amb_min           = W_AMBIENT_MIN,
  w_amb_max           = W_AMBIENT_MAX,
  w_amb_p             = W_AMBIENT_P,
  reef_radius_m       = REEF_RADIUS_M,
  reef_radius_px      = REEF_RADIUS_PX_TARGET,
  maxsearch_m         = MAXSEARCH_M,
  ambient_maxsearch_m = AMBIENT_MAXSEARCH_M,
  ambient_sunbreaks   = AMBIENT_SUNBREAKS,
  timestep            = TIMESTEP_INTERVAL,
  buffer_m            = BUFFER_M,
  pool_sites          = paste(pool_keys, collapse = ";"),
  n_hourly            = hourly_overall$n,
  rmse_hourly         = hourly_overall$rmse,
  r2_hourly           = hourly_overall$r2,
  r2_1to1_hourly      = hourly_overall$r2_1to1,
  bias_hourly         = hourly_overall$bias,
  n_dli               = dli_overall$n,
  rmse_dli            = dli_overall$rmse,
  r2_dli              = dli_overall$r2,
  r2_1to1_dli         = dli_overall$r2_1to1,
  bias_dli            = dli_overall$bias,
  stringsAsFactors    = FALSE
)

log_path <- file.path(RUNS_ROOT, "metrics_log.csv")
if (file.exists(log_path)) {
  write.table(log_row, log_path, sep = ",", append = TRUE,
              row.names = FALSE, col.names = FALSE)
} else {
  write.csv(log_row, log_path, row.names = FALSE)
}
message("Appended sweep summary to ", log_path)

# ── save per-run validation CSV ───────────────────────────────────────────────
val_path <- file.path(RUNS_DIR, sprintf("validation_%s.csv", RUN_TAG))
write.csv(validation, val_path, row.names = FALSE)
message("Saved ", val_path)


# -- depth bias ---------------------------------------------------------------

val_with_depth <- validation %>%
  mutate(residual = pred_par - light_mol,
         depth_bin = cut(water_depth_m, 
                         breaks = c(0, 2, 2.5, 3, 3.5, 4, 5),
                         labels = c("<2m","2-2.5","2.5-3","3-3.5","3.5-4",">4m")))

bias_by_depth <- val_with_depth %>%
  group_by(sitekey, depth_bin) %>%
  summarise(n = n(),
            mean_meas = mean(light_mol, na.rm=TRUE),
            mean_pred = mean(pred_par, na.rm=TRUE),
            bias_abs  = mean(residual, na.rm=TRUE),
            bias_rel  = mean(residual/light_mol, na.rm=TRUE) * 100,
            rmse      = sqrt(mean(residual^2, na.rm=TRUE)),
            .groups = "drop")

print(bias_by_depth, n=Inf)

# ══════════════════════════════════════════════════════════════════════════════
#  OPTIONAL — FIT W_DIRECT / W_AMBIENT (and Kd). Pools ALL site/years in
#  `validation`, so multi-site just gives it a bigger fitting set.
#
#  CAVEAT in altitude-dependent mode: this searches a STATIC scalar, but the
#  pipeline is using a curve. The "best static" reported here is informational
#  — it's the best constant that would approximate your curve on this data —
#  NOT a recommendation to switch back to static mode. Warning prints to remind.
# ══════════════════════════════════════════════════════════════════════════════
fit_ambient_weights <- function(validation,
                                kd_grid = seq(0.10, 0.60, by = 0.02),
                                wamb_grid = seq(0, 0.95, by = 0.05)) {
  if (USE_ALT_DEPENDENT_WEIGHTS) {
    warning("fit_ambient_weights() searches a STATIC W_AMBIENT, but ",
            "USE_ALT_DEPENDENT_WEIGHTS=TRUE. The 'best' value here is the ",
            "constant that best approximates the altitude curve — it does NOT ",
            "tell you to disable the altitude curve. To genuinely co-fit the ",
            "altitude endpoints, expand to a 3-D grid (Kd x W_AMBIENT_MIN x ",
            "W_AMBIENT_MAX) and apply w_ambient_alt() inside the loop.")
  }
  fitdat <- validation %>%
    filter(!is.na(shadow_fraction), !is.na(ambient_fraction),
           !is.na(light_mol), !is.na(surface_par), !is.na(water_depth_m))
  
  rmse_grid <- matrix(NA_real_, length(kd_grid), length(wamb_grid),
                      dimnames = list(paste0("kd", kd_grid),
                                      paste0("wamb", wamb_grid)))
  for (a in seq_along(kd_grid)) for (b in seq_along(wamb_grid)) {
    wa  <- if (USE_AMBIENT) wamb_grid[b] else 0
    eff <- (1 - wa) * fitdat$shadow_fraction + wa * fitdat$ambient_fraction
    pred <- fitdat$surface_par * eff * exp(-kd_grid[a] * fitdat$water_depth_m)
    rmse_grid[a, b] <- sqrt(mean((pred - fitdat$light_mol)^2, na.rm = TRUE))
  }
  
  bi <- which(rmse_grid == min(rmse_grid), arr.ind = TRUE)[1, ]
  best_kd   <- kd_grid[bi[1]]
  best_wamb <- wamb_grid[bi[2]]
  
  cat("\n=====================================================\n")
  cat("AMBIENT-WEIGHT CO-FIT (geometry-based diffuse layer)\n")
  cat("=====================================================\n")
  cat("Best Kd            :", best_kd, "\n")
  cat("Best W_AMBIENT     :", best_wamb, " (W_DIRECT =", 1 - best_wamb, ")\n")
  cat("Best RMSE          :", round(min(rmse_grid), 1), "\n")
  
  dev.new(width = 8, height = 6)
  par(mar = c(5, 5, 4, 6))
  image(kd_grid, wamb_grid, rmse_grid, col = hcl.colors(100, "viridis"),
        xlab = "Kd (m^-1)", ylab = "W_AMBIENT (diffuse weight)",
        main = "RMSE: Kd vs ambient weight")
  contour(kd_grid, wamb_grid, rmse_grid, add = TRUE, col = "white", nlevels = 12)
  points(best_kd, best_wamb, pch = 4, col = "red", cex = 2, lwd = 3)
  
  invisible(list(best_kd = best_kd, best_wamb = best_wamb,
                 rmse = min(rmse_grid), grid = rmse_grid))
}


fit <- fit_ambient_weights(validation)

# usage after pipeline (only meaningful if USE_ALT_DEPENDENT_WEIGHTS = FALSE):
#   fit <- fit_ambient_weights(validation)
#   W_DIRECT  <- 1 - fit$best_wamb
#   W_AMBIENT <- fit$best_wamb
#   # then re-run run_pipeline() to apply the fitted weights