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
#  capture_tide_m is the tide at each site's photogrammetry timestamp, linearly
#  interpolated from the SAME env_long$tide_cm series the predictions use, so the
#  correction is internally consistent (one tide reference throughout).
#
#  NOTE ON KEYS: env_long / geom use site = "TB"/"CB"/"RS" and year = 18/19, so
#  the photogrammetry table below is keyed the SAME way. If it were keyed on the
#  long names ("turtle_beach", 2018) the left_join would silently yield NA and
#  every water_depth_m would become NA.
# ══════════════════════════════════════════════════════════════════════════════

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
#  REPLACES the old approx()/env_long clamp. That version snapped to the nearest
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
  
  # extrapolation check
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
  # "CB18" = 0.00   # <-- REQUIRED for your run: external Lizard tide @ 2018-11-02 08:32
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


# ── FILE-PATH ROOTS ───────────────────────────────────────────────────────────
#create filepath for DEMS
DEM_DIR <- "/Users/samueleszenyi/Desktop/MEM MSc/1thesis/dems"

# ── constants ─────────────────────────────────────────────────────────────────

KD_CONSTANT       <- 0.2  #fallback value      # KD Constant — universal diffuse attenuation coefficient (m⁻¹)

# which KdPAR variant to apply:
#   "kdpar_photon" = photon-weighted full-band integral 
#   "kdpar_mean"   = unweighted mean of the 10 bands    
#   "kdpar_geom"   = geometric mean, red-outlier-resistant
KDPAR_COL <- "kdpar_photon"

kdpar_lookup <- read.csv("kdpar_lookup.csv", stringsAsFactors = FALSE)
# build a named vector keyed "SITEyy" -> chosen KdPAR, e.g. "CB18" -> 0.2358
kdpar_lookup$key <- paste0(kdpar_lookup$site, sprintf("%02d", kdpar_lookup$year %% 100))
KDPAR_BY_KEY <- setNames(kdpar_lookup[[KDPAR_COL]], kdpar_lookup$key)

message("Loaded per-site KdPAR (", KDPAR_COL, ") for ", length(KDPAR_BY_KEY),
        " site/years. Fallback Kd = ", KD_CONSTANT)

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
# TRUE  = derive Kd(t) per site by inverting Beer-Lambert on ONE open, known-depth
#         reference logger: Kd(t) = -ln(ref_light / surface_par) / ref_depth.
#         The reference logger is an INPUT (like surface_par) and is EXCLUDED
#         from validation, so this is not circular — but it requires a genuinely
#         UNSHADED reference at correct depth, and ref_light & surface_par in the
#         SAME units. CHECK the printed transmittance diagnostics before trusting.
USE_REFERENCE_KD <- TRUE

# One reference hoboID per sitekey. MUST be an open/unshaded logger.
# Pick these from your shadow maps — the logger with the highest, flattest shadow
# (near 1) and no overhang. Leave a site out to fall back to satellite Kd for it.
REFERENCE_LOGGER <- c(
  "TB18" = "H34"
)



TIMESTEP_INTERVAL <- "1 hour"    # the timestep interval to measure

# ── DOWNSAMPLING: target output size, not a fixed factor ──────────────────────
# OLD behaviour: a fixed 10x factor. Problem — output resolution then depended on
# native DEM size, so a 3000px DEM became 300px (too coarse, reef structure lost)
# while a 15000px DEM became 1500px. The 300px ones badly degrade the shadows.
#
# NEW behaviour: pick a TARGET pixel width and derive the factor PER DEM so every
# DEM ends up near the same working resolution. Factor is computed in prepare_dem
# as round(native_px / TARGET), floored at 1 (never UPsample). A DEM already at or
# below the target is left at its native resolution (factor 1).
TARGET_PIXELS     <- 1000    # desired output width/height in pixels after downsample
MIN_DOWNSAMPLE    <- 1       # never upsample; 1 = keep native if already small
# zscale still adapts automatically — it is pixel_size * (the per-DEM factor),
# so coarser or finer DEMs each get the correct vertical scaling.
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

#Buffer to extract shadow and ambient values (always matchign) in meters
BUFFER_M <- 0.01

MAXSEARCH_FACT    <- 700    # change maxsearch distance, default 200
# every pixel shoots a ray up towards the sun at x distance, this is in pixel units, so if the DEM diameter is 13m, current pixel size = 0.02m, 
#"the maximum horizontal distance that the system should propagate rays." - rayshader documentation 
# if it hits something within that distance, it is a shadow

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

AMBIENT_SUNBREAKS  <- 50    # rays sent in a circle around each pixel (rayshader default is 24)
AMBIENT_MAXSEARCH  <- 500    # how far (in PIXELS) to look for sky-blocking geometry.
# default is 30; on a downsampled reef DEM 50–100 is
# reasonable; larger = sees more distant blockers = slower.
W_DIRECT           <- 0.4  # weight on the sun-beam (ray_shade) component
W_AMBIENT          <- 0.6  # weight on the diffuse (ambient_shade) component
# NB: W_DIRECT + W_AMBIENT must equal 1
# (fit these with fit_ambient_weights() at the end)

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
  if (USE_REFERENCE_KD)   return("reference Kd")
  if (USE_SPECTRAL_KD)    return("spectral Kd")
  paste0("satellite Kd (", KDPAR_COL, ")")
}
KD_METHOD <- kd_method_label()   # evaluated once; re-call kd_method_label() if you flip toggles
print(KD_METHOD)

# ══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# reusable code chunk recipes that I made to be called later
# ══════════════════════════════════════════════════════════════════════════════

# ── prepare_dem ───────────────────────────────────────────────────────────────
# Load DEM, reproject to UTM, downsample, convert to a matrix for rayshading.
# DEMs are originally in WGS84 geographic degrees (EPSG:4326).
# Reprojecting to UTM 55S puts x, y, z all in metres so zscale = 1.
# NOW ALSO computes the ambient occlusion matrix once and ships it in the list.
prepare_dem <- function(site, year, target_px = TARGET_PIXELS) {
  path <- file.path(DEM_DIR, paste0(site, year, ".tif"))
  message("  Loading DEM: ", path)
  dem  <- rast(path) |> trim()
  message("  Raw DEM CRS: ", crs(dem, describe=TRUE)$name)
  
  # reproject to UTM zone 55S so x, y, z are all in metres
  # this means zscale = 1 (1 horizontal metre = 1 vertical metre)
  dem  <- project(dem, UTM_CRS)
  message("  Reprojected to UTM 55S (EPSG:32755)")
  message("  DEM extent (metres): ", paste(round(as.vector(ext(dem))), collapse = " "))
  
  # ── per-DEM downsample factor: aim for ~target_px on the larger dimension ────
  # aggregate(fact = N) divides each side by N, so to land near target_px we use
  # native_px / target_px. floor-guarded at MIN_DOWNSAMPLE so we never upsample,
  # and a DEM already <= target keeps factor 1 (native resolution).
  native_px <- max(nrow(dem), ncol(dem))
  fact      <- max(MIN_DOWNSAMPLE, round(native_px / target_px))
  message("  Native size: ", nrow(dem), " x ", ncol(dem),
          " (max ", native_px, " px) -> downsample factor ", fact,
          " -> ~", round(native_px / fact), " px target")
  
  dem_small  <- aggregate(dem, fact = fact)
  
  # zscale = xy pixel spacing / z (all UTM metres, so xy = z)
  res_x      <- (ext(dem)$xmax - ext(dem)$xmin) / ncol(dem)
  res_y      <- (ext(dem)$ymax - ext(dem)$ymin) / nrow(dem)
  pixel_size <- (res_x + res_y) / 2
  zscale     <- pixel_size * fact
  message("  Pixel size (UTM metres): ", round(pixel_size, 6))
  message("  Zscale (downsampled): ",    round(zscale, 6))
  
  # ── ALL masking/filling done on the RASTER, so dem_small and the matrix
  #    derived from it share one identical footprint, extent and orientation.
  #
  #  WHY THIS CHANGED
  #  ────────────────
  #  The old mask drew a circle of radius 0.90*min(nr,nc)/2 centred on the
  #  BOUNDING-BOX centre. But the real survey (a ~13 m blob) is NOT centred in
  #  the bounding box, and the off-reef area after reprojection is a ROTATED
  #  SQUARE of real-valued fill, not NA. So the circle:
  #    • bled into the rotated-square fill on some sides  -> dark corner wedges
  #      (the CB19/year-19 "unclipped" look: rotated square showing through), and
  #    • clipped real reef on the narrow side of off-centre / rectangular DEMs
  #      (the GT19 rectangular problem).
  #
  #  THE FIX (geometry-only, works for every site without tuning)
  #  ────────────────────────────────────────────────────────────
  #  1. Find where the DEM ACTUALLY has data (data_present = !is.na).
  #  2. Centre the circle on the DATA centroid, not the bbox centre.
  #  3. Size the radius from the data's own spread (a high quantile of pixel
  #     distances from that centroid), so it rounds off ragged edges into a
  #     clean disc that fits inside the real survey.
  #  4. INTERSECT the circle with data_present. This is the key line: a pixel is
  #     only "inside" if it is BOTH within the disc AND real data, so the
  #     rotated-square fill can never enter the reef region again.
  #  All four DEMs (clipped or unclipped, rotated square or rectangle) collapse
  #  to the same clean circular reef.
  
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
  
  # radius from the data's OWN spread: 95th percentile of real-pixel distances,
  # so a few stray edge pixels don't inflate it. The disc then sits just inside
  # the real survey rather than reaching into the bounding-box corners.
  dist_all <- sqrt((rr - row_c)^2 + (cc - col_c)^2)
  rad      <- as.numeric(stats::quantile(dist_all[data_present], 0.95)) * 1.00
  
  circle_mask <- dist_all <= rad
  
  # KEY INTERSECTION: inside == within the disc AND genuinely real data.
  inside_mat  <- circle_mask & data_present
  message("  Reef mask: ", sum(inside_mat), " px inside ",
          "(centroid row/col ", round(row_c), "/", round(col_c),
          ", radius ", round(rad), " px)")
  
  # write back: everything outside the intersected mask -> NA (row-major vector)
  v <- v0
  v[!as.vector(t(inside_mat))] <- NA
  values(dem_small) <- v
  
  # fill the off-reef area BELOW the real reef minimum (not equal to it),
  # so the border sits well under everything and cannot cast shadows inward.
  reef_min  <- min(values(dem_small), na.rm = TRUE)
  reef_max  <- max(values(dem_small), na.rm = TRUE)
  fill_val  <- reef_min - (reef_max - reef_min)   # one full relief-depth below the reef
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
       ambient   = amb_matrix)         # <-- NEW: diffuse occlusion layer
}
#this returns: a matrix for rayshader, a raster for coordinate extraction,
#the zscale, the reef mask, and now the ambient occlusion matrix.

# ── make_ambient ──────────────────────────────────────────────────────────────
#  Compute the ambient occlusion matrix ONCE per DEM. Same orientation/units as
#  ray_shade so it samples through the exact same extract path.
#  zscale must match what you pass to ray_shade (= dem$zscale, UTM metres).
make_ambient <- function(dem_matrix, zscale,
                         sunbreaks = AMBIENT_SUNBREAKS,
                         maxsearch = AMBIENT_MAXSEARCH) {
  message("  Computing ambient occlusion (sunbreaks=", sunbreaks,
          ", maxsearch=", maxsearch, ") ...")
  amb <- ambient_shade(
    heightmap = dem_matrix,
    sunbreaks = sunbreaks,
    maxsearch = maxsearch,
    zscale    = zscale,
    multicore = FALSE # set TRUE to use all cores; needs the parallel backend
  )
  message("    ambient range: ", round(min(amb, na.rm = TRUE), 3),
          " to ",               round(max(amb, na.rm = TRUE), 3),
          " | mean: ",          round(mean(amb, na.rm = TRUE), 3))
  amb
}

# ── load_hobo_locations ───────────────────────────────────────────────────────
# build a sf point object from geom_new coordinates for a given site/year.
# coords are in the local transverse mercator projection used when the
# geom script vivi made extracted them from the DEMs.

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
#0 is fully shaded, and 1 is fully lit, if the sun is below the horizon it skips the measurements to save time
#zscale = 1 because DEM is in UTM metres (x, y, z all same units)
make_shadow <- function(dem_matrix, sun, zscale) {
  if (sun$altitude <= 0.05) {
    message("    (sun below horizon — skipping raytrace)")
    return(NULL)
  }
  
  shadow <- ray_shade(
    dem_matrix,
    sunaltitude = sun$altitude * 180 / pi, #gets sun altitude in degrees from radians
    sunangle    = sun$azimuth  * 180 / pi, #gets sun azimuth in degrees from radians
    zscale      = zscale,
    maxsearch   = MAXSEARCH_FACT #limits how far each ray travels
  )
  
  message("    shadow range: ", round(min(shadow, na.rm = TRUE), 3),
          " to ",               round(max(shadow, na.rm = TRUE), 3),
          " | mean: ",          round(mean(shadow, na.rm = TRUE), 3))
  shadow
}


# ── calc_water_attenuation ────────────────────────────────────────────────────
# Beer-Lambert law: exp(−Kd × depth). Single-KdPAR path (default).
# Returns the fraction of surface PAR surviving to depth; deeper logger = smaller.
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
  if (is.na(water_depth_m) || water_depth_m <= 0) return(NA_real_)   # <- add is.na() guard
  ok <- is.finite(kd_bands)
  if (sum(ok) < 2) return(NA_real_)
  kb <- kd_bands[ok]; wv <- waves[ok]
  f  <- wv / sum(wv)
  sum(f * exp(-kb * water_depth_m))
}


# ── derive_reference_kd ───────────────────────────────────────────────────────
# For one site/year, invert Beer-Lambert on the reference logger to get Kd per
# hour. Returns a named vector keyed by the hourly POSIXct (as character) -> Kd.
# Kd(t) = -ln( ref_light(t) / surface_par(t) ) / ref_depth(t)
# Guards: needs ref_light < surface_par (transmittance in (0,1)), depth > 0.
derive_reference_kd <- function(env_site, ref_id, cap_tide_m) {
  ref <- env_site %>%
    filter(as.character(hoboID) == ref_id) %>%
    mutate(datetime_pos = floor_date(datetime_pos, TIMESTEP_INTERVAL)) %>%
    group_by(datetime_pos) %>%
    summarise(surface_par = mean(surface_par, na.rm = TRUE),
              ref_light   = mean(light_mol,   na.rm = TRUE),
              tide_cm     = mean(tide_cm,     na.rm = TRUE),
              dpHB        = mean(dpHB,         na.rm = TRUE), .groups = "drop") %>%
    mutate(ref_depth = (tide_cm / 100) - dpHB - cap_tide_m,
           transmit  = ref_light / surface_par,
           Kd_t      = ifelse(transmit > 0 & transmit < 1 & ref_depth > 0,
                              -log(transmit) / ref_depth, NA_real_))
  
  # DIAGNOSTIC — you MUST eyeball this. If transmit is >1 or absurd, units differ.
  message("  [ref-Kd] ", ref_id, " transmit range ",
          round(min(ref$transmit, na.rm=TRUE), 3), "–",
          round(max(ref$transmit, na.rm=TRUE), 3),
          " | Kd(t) range ",
          round(min(ref$Kd_t, na.rm=TRUE), 3), "–",
          round(max(ref$Kd_t, na.rm=TRUE), 3),
          " | median Kd ", round(median(ref$Kd_t, na.rm=TRUE), 3))
  
  setNames(ref$Kd_t, as.character(ref$datetime_pos))
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
  
  # buffer in meters dist = 0.03 = 0.03m (3cm)
  hobo_buffered <- st_buffer(hobo_utm, dist = BUFFER_M)
  
  # extract mean shadow value within each buffer
  vals <- terra::extract(shadow_rast, vect(hobo_buffered), fun = mean, na.rm = TRUE)
  vals[, 2]
}

# ── extract_hobo_ambient ──────────────────────────────────────────────────────
#  Sample the ambient occlusion matrix at each logger, using identical machinery
#  to extract_hobo_shadow so the two layers line up pixel-for-pixel.
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
#the core calculations for each timestep (I chose 1 hour)
#  NOW blends a directional (ray_shade) and a diffuse (ambient_shade) layer.
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
  # zscale = 1 because DEM is in UTM metres
  shadow <- make_shadow(dem_matrix, sun, zscale = dem_zscale)
  # blank off-reef border so it can't be sampled (guard: shadow is NULL at night)
  if (!is.null(shadow)) shadow[!dem_reef_mask] <- NA
  
  # 2b. ambient occlusion (diffuse) — sun-independent, blank border the same way
  ambient_masked <- dem_ambient
  ambient_masked[!dem_reef_mask] <- NA
  
  # 3. sample BOTH layers at each hobo logger.
  #if the sun is below the horizon (dark hours) shadow is NA and skipped;
  #ambient is geometry-only so it is always sampled.
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
  message("    dem_rast extent: ", paste(round(as.vector(ext(dem_rast))), collapse = " "))
  hobo_proj <- st_transform(hobo_sf, crs(dem_rast))
  message("    hobo bbox: ", paste(round(st_bbox(hobo_proj)), collapse = " "))
  message("    dem_rast CRS: ", crs(dem_rast, describe=TRUE)$name)
  message("    hobo_sf CRS: ", st_crs(hobo_sf)$input)
  
  # 4. this makes a table pairing each shadow + ambient value with the logger id
  shadow_df <- data.frame(
    hoboID_env       = hobo_sf$hoboID,
    shadow_fraction  = shadow_vals,
    ambient_fraction = ambient_vals
  )
  
  # 5. this joins the fractions with logger data (depth, surface par, kd), blends
  # the directional + diffuse layers, then computes the final predicted value as
  # surface par x (w_direct·shadow + w_ambient·ambient) x water attenuation
  result <- shadow_df %>%
    left_join(
      env_block %>%
        mutate(hoboID_env = as.character(hoboID)) %>%
        select(hoboID_env, Kd, water_depth_m, surface_par, site, year),
      by = "hoboID_env"
    ) %>%
    mutate(
      datetime_pos = timestamp,
      sun_altitude = sun$altitude,
      sun_azimuth  = sun$azimuth,
      # water attenuation: single-Kd (Beer-Lambert) OR per-band spectral sum.
      # rowwise over water_depth_m because spectral atten is depth-specific.
      water_atten  = if (USE_SPECTRAL_KD && !is.null(kd_bands)) {
        vapply(water_depth_m,
               function(z) calc_water_attenuation_spectral(kd_bands, WAVES_PAR, z),
               numeric(1))
      } else {
        calc_water_attenuation(Kd, water_depth_m)
      },
      eff_light    = if (USE_AMBIENT) {
        W_DIRECT * shadow_fraction + W_AMBIENT * ambient_fraction
      } else {
        shadow_fraction            # ambient off: direct beam only
      },      pred_par     = surface_par * eff_light * water_atten
    )
  
  as_tibble(result)
}


# ══════════════════════════════════════════════════════════════════════════════
#  main code loop that puts it all together
# master function that does following:

# loops over every site x year combination in data
# loads elevation map and logger positions ONCE
#builds a schedule of hourly daytime timesteps to process
#and for every hour that qualifies (sun over horizon), raytraces DEM and computes PAR at every logger location

#the outer loop sets up the machinery for every site/year, and the inner loop runs it one hour at a time

#inputs, env_long - full environmental data table
#geom - logger coordinate lookup table
#dems

#output - one table with one row per logger per timestep, with pred_par
# ══════════════════════════════════════════════════════════════════════════════

run_pipeline <- function(env_long, geom, dem_store = NULL) {
  #results in an empty list that it fills one bit at a time
  #dem_store (optional): named list keyed "TB18" etc, built ONCE above. If
  #supplied, the spatial data (DEM, ambient occlusion, hobo points) is REUSED
  #from it instead of calling prepare_dem() again — prepare_dem (esp.
  #ambient_shade) is the slow step and was previously run twice per site/year.
  
  results <- list()
  #get the unique site codes present (CB, L1, Gt, etc)
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
    #when they are outside hourly loop we only do this once, since it takes a lot of time (i broke my R doing it inside the loop)
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
      
      # ── capture-tide for THIS site/year (replaces the old hardcoded "- 2") ──
      # one value, looked up from the capture_tide table built near the top.
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
      ref_kd_lookup <- NULL
      if (USE_REFERENCE_KD) {
        ref_id <- REFERENCE_LOGGER[[make_key(s, y)]]
        if (!is.null(ref_id)) {
          env_site_full <- env_long %>% filter(site == s, year == y)
          ref_kd_lookup <- derive_reference_kd(env_site_full, ref_id, cap_tide_m)
          message("  Using reference logger ", ref_id, " for time-varying Kd")
        } else {
          message("  No reference logger for ", make_key(s, y),
                  " — falling back to satellite Kd")
        }
      }
      
      #---------step2 build hourly schedule of timesteps
      #the raw data is 5 minute intervals, but we only want to raytrace 1 hour to save computing
      #the section here collapses 5 min data to one row per logger per hour, computes sun angle, and throws away hours where the sun is too low to mean anything
      #%>% chains everything together in order
      hourly_times <- env_long %>%
        
        #keep only rows belonging to this site year
        filter(site == s, year == y) %>%
        
        #sort chronologically for floor date to work
        arrange(datetime_pos) %>%
        
        #snap every 5 min timestamp to its hour
        #11:05, 11:10, 11:55 all become 11, which makes rows share the same datetime
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
        
        #give every row our constant KD value, the kd values from env_long are discarded, since using them would be circular
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
          water_depth_m = (tide_cm / 100) - dpHB - cap_tide_m
        ) %>%
        
        # Each site/year now carries its own KdPAR. The message lines that print
        # KD_CONSTANT will still print the fallback value, which is fine — the actual Kd
        # used is in env_block$Kd and already prints via process_timestep's
        # "| Kd: round(env_block$Kd[1], 4)" message, so you can confirm per-site.
        
        
        # tide-corrected water column above sensor (recomputed here, not upstream)
        # tide_cm/100 = tide above datum (m); dpHB = sensor bed elevation (m, negative down)
        # cap_tide_m  = actual tide at this site's photogrammetry capture time,
        #               interpolated from env_long$tide_cm (replaces the old flat -2).
        
        #r applies the next mutate one row at a time, getsunlightposition needs one single timestamp
        rowwise() %>%
        
        #compute solar altitude in degrees, stores it as sun_alt, converting radians to degrees by altitude*180/pi
        #centroid is in WGS84 lat/lon as required by suncalc
        mutate(
          sun_alt = getSunlightPosition(
            datetime_pos, centroid$lat, centroid$lon
          )$altitude * 180 / pi #return degrees
        ) %>%
        ungroup() %>%
        
        #discard hours where the sun is less than 15 degrees over the horizon
        #at super low angles, the shadows hit dem imperfections harder and the raytracer produces unreliable shadows
        #we lose a bit of early morning data but it can always be changed
        filter(sun_alt > 20) %>%
        
        #remove sun_alt column, since it served its purpose
        select(-sun_alt)
      
      # DIAGNOSTIC — confirm depth varies with tide; remove after checking
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
      #for each hour i pull out all logger rows for that hour, pass them to process_timestep, which raytraces once then extracts shadow values, then output results
      
      #for seq_len(n_steps) produces 1,2,3,n_steps, which is safer than writing 1:n_steps because if n_steps were 0, 1:0 would loop twice
      for (i in seq_len(n_steps)) {
        
        #pick the i-th timestamp from the hours we decided
        t         <- timestamps[i]
        
        #pull out logger rows for this hour, [] is a row filter, which keeps rows where datetime_pos == t
        #process_timestep needs all loggers at once so it can extract shadow values for each one from the shadow map that was made earlier
        env_block <- hourly_times[hourly_times$datetime_pos == t, ]
        
        message("\n[", i, "/", n_steps, "] site: ", s, " | year: ", y,
                " | loggers in block: ", nrow(env_block),
                " | Kd: ", round(env_block$Kd[1], 4))   # was KD_CONSTANT
        
        #run the full prediction for the timestep
        #solar position -> raytrace + ambient -> logger shadows -> beer lambert law -> predicted par (pred_par)
        # per-site 10-band Kd vector for the spectral path (NULL-safe if absent)
        site_kd_bands <- KD_BANDS_BY_KEY[[paste0(s, sprintf("%02d", as.integer(y) %% 100))]]
        
        res <- process_timestep(
          env_block,   #data frame: one row per logger for this hour
          dem$matrix,  #dem as a numeric matrix
          dem$rast,    #dem as a spatial object/raster for coordinate stuff (UTM)
          dem$zscale,  #zscale = 1, x/y/z all in metres after UTM reprojection
          hobo_sf,     #logger locations as a spatial object
          centroid,    #site lat/lon (WGS84) for suncalc sun angle lookup
          dem$reef_mask,  # blank the off-reef border before sampling
          dem$ambient,    # precomputed ambient occlusion matrix
          site_kd_bands   # <-- NEW: per-band Kd for spectral attenuation
        )
        
        #add hour's result tibble to results
        #length(results) +1 puts it as the next empty slot
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
#  HOW MULTI-SITE WORKS NOW
#  ────────────────────────
#  Nothing below hardcodes a single site. The flow is:
#    1. `wanted`     : the site/year combos you want this run.
#    2. run_pipeline : loops them internally, returns one pooled `output`.
#    3. dem_store    : DEMs/hobo points are loaded ONCE here and cached in a
#                      named list keyed "TB18", "CB19", ... so every plot reuses
#                      them instead of calling prepare_dem again.
#    4. per-site figs: PAR-per-logger panels and shadow maps LOOP over each
#                      site/year and save one file each.
#    5. pooled figs  : DLI + time-of-day scatter pool all rows; point SHAPE = the
#                      site/year, COLOUR = depth, and every sensor is labelled.
#
#  The one rule that makes it all work: a logger is only unique WITHIN a
#  site/year, so everywhere we used hoboID alone we now use a sitekey + hoboID.
# ══════════════════════════════════════════════════════════════════════════════

# ── helper: one canonical site/year key, used everywhere ──────────────────────
make_key <- function(site, year) paste0(site, year)   # "TB" + 18 -> "TB18"

# ── FILTER DATA ───────────────────────────────────────────────────────────────
#  Edit ONLY this table. One row per site/year you want in the run.
#  (A single-site run is just a one-row table — no other code changes needed.)
wanted <- tibble::tribble(
  ~site, ~year,
  "TB", 18L
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

# the distinct site/year combos actually present after filtering (data-driven)
site_years <- env_long_test %>%
  distinct(site, year) %>%
  mutate(sitekey = make_key(site, year)) %>%
  arrange(sitekey)

message("Run covers ", nrow(site_years), " site/year combos: ",
        paste(site_years$sitekey, collapse = ", "))


# ══════════════════════════════════════════════════════════════════════════════
#  LOAD ALL SPATIAL DATA ONCE  →  dem_store
#  Replaces both the old single-site "DEM VERIFICATION" and "save to GE" blocks.
#  Each entry holds everything the plots need for that site/year.
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
#  Set to FALSE to skip once you trust the DEMs.
VERIFY_DEMS <- FALSE
if (VERIFY_DEMS) {
  for (key in names(dem_store)) {
    e  <- dem_store[[key]]
    dm <- e$dem
    
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
# Drop reference loggers from validation when reference-Kd is on (they're inputs).
if (USE_REFERENCE_KD) {
  ref_drop <- tibble::tibble(
    sitekey    = names(REFERENCE_LOGGER),
    hoboID_env = as.character(REFERENCE_LOGGER))
  validation <- validation %>% anti_join(ref_drop, by = c("sitekey", "hoboID_env"))
  message("Validation rows after dropping reference loggers: ", nrow(validation))
}

message("Validation rows: ", nrow(validation),
        " across ", length(unique(validation$sitekey)), " site/years")

# ══════════════════════════════════════════════════════════════════════════════
#  REFERENCE-LOGGER FINDER  (constant, output-driven recommendation)
#  ──────────────────────────────────────────────────────────────────────────
#  Picks ONE representative logger per site to invert Beer-Lambert for Kd(t).
#  The ideal reference is, across the whole daytime cycle:
#    • as UNSHADED as possible            -> high mean shadow_fraction
#    • as STABLY unshaded as possible     -> high MINIMUM shadow + low variance
#                                            ("open throughout the day", not just
#                                             bright at noon and shaded at the
#                                             ends — those give a biased Kd)
#    • as OPEN TO THE SKY as possible      -> high mean ambient_fraction (least
#                                            ambient occlusion, sees full sky)
#    • a VALID, STABLE known depth         -> hard gate (needed to invert B-L)
#    • ENOUGH daytime hours to trust it    -> coverage gate (n_hours)
#
#  WHY SCORE FROM `output`, NOT `validation`
#  ─────────────────────────────────────────
#  validation DROPS the reference loggers when USE_REFERENCE_KD = TRUE, so
#  scoring from validation is circular: flip the toggle / re-run and the
#  candidate pool changes, so the recommended logger can change. `output` is the
#  full pre-validation prediction table (every logger, every hour, never
#  dropped) and depends ONLY on geometry + depth, NOT on the Kd method or on
#  which logger was previously chosen. Scoring from it gives a STABLE, constant
#  recommendation you can paste once into REFERENCE_LOGGER and leave fixed.
#  (It also doesn't need light_mol, so it's independent of the validation join.)
#
#  NB: still a RECOMMENDATION — verify the top pick against the shadow maps
#  (is it actually on open sand, not perched on a coral head?).
# ══════════════════════════════════════════════════════════════════════════════
find_reference_loggers <- function(out,
                                   depth_sd_max = 3,   # max hour-to-hour depth sd (m)
                                   min_hours    = 5,   # min daytime hours to qualify
                                   w_shadow_mean = 0.5,  # weight: avg openness to sun
                                   w_shadow_min  = 0.25,  # weight: worst-case openness
                                   w_shadow_stab = 0.15,  # weight: temporal steadiness
                                   w_ambient     = 0.1) {# weight: openness to sky
  stopifnot(abs(w_shadow_mean + w_shadow_min + w_shadow_stab + w_ambient - 1) < 1e-9)
  
  scored <- out %>%
    filter(!is.na(shadow_fraction)) %>%      # daytime, raytraced rows only
    group_by(sitekey, hoboID_env) %>%
    summarise(
      shadow_mean  = mean(shadow_fraction,  na.rm = TRUE),
      shadow_min   = min(shadow_fraction,   na.rm = TRUE),  # worst hour of the day
      shadow_sd    = sd(shadow_fraction,    na.rm = TRUE),  # day-cycle variability
      ambient_mean = mean(ambient_fraction, na.rm = TRUE),
      depth_mean   = mean(water_depth_m,    na.rm = TRUE),
      depth_sd     = sd(water_depth_m,      na.rm = TRUE),
      n_hours      = dplyr::n_distinct(datetime_pos),
      .groups = "drop"
    ) %>%
    # HARD GATES: usable depth to invert Beer-Lambert, and enough daytime cover.
    mutate(
      depth_ok = is.finite(depth_mean) & depth_mean > 0 &
        is.finite(depth_sd)   & depth_sd < depth_sd_max,
      hours_ok = n_hours >= min_hours,
      gate_ok  = depth_ok & hours_ok
    )
  
  # Normalise each openness axis to 0–1 WITHIN site so weights mean the same
  # thing everywhere (each site's own best logger sits at ~1 per axis).
  # shadow_sd is inverted (low variance = steady = good -> high stability_norm).
  norm01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    out <- (x - r[1]) / (diff(r) + 1e-9)
    out[!is.finite(out)] <- 0   # non-finite (e.g. sd of 1 obs) -> worst, never NA
    out
  }
  scored <- scored %>%
    group_by(sitekey) %>%
    mutate(
      shadow_mean_norm = norm01(shadow_mean),
      shadow_min_norm  = norm01(shadow_min),
      stability_norm   = 1 - norm01(shadow_sd),   # invert: steadier -> higher
      ambient_norm     = norm01(ambient_mean),
      open_score = ifelse(
        gate_ok,
        w_shadow_mean * shadow_mean_norm +
          w_shadow_min  * shadow_min_norm  +
          w_shadow_stab * stability_norm   +
          w_ambient     * ambient_norm,
        0)
    ) %>%
    arrange(sitekey, desc(open_score)) %>%
    ungroup()
  
  # ranked table per site
  for (k in sort(unique(scored$sitekey))) {
    cat("\n── Reference candidates for ", k, " (ranked) ──\n", sep = "")
    sub <- scored %>% filter(sitekey == k) %>%
      transmute(hoboID_env,
                shadow      = round(shadow_mean, 3),
                shadow_min  = round(shadow_min,  3),
                shadow_sd   = round(shadow_sd,   3),
                ambient     = round(ambient_mean, 3),
                depth_m     = round(depth_mean,  2),
                depth_sd    = round(depth_sd,    2),
                n_hours,
                gate_ok,
                score       = round(open_score,  3))
    print(as.data.frame(sub), row.names = FALSE)
  }
  
  # best per site (highest score, must pass both gates)
  best <- scored %>%
    filter(open_score > 0) %>%
    group_by(sitekey) %>%
    slice_max(open_score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(sitekey, hoboID_env, shadow_mean, shadow_min,
           ambient_mean, depth_mean, n_hours, open_score)
  
  # flag any site that produced NO qualifying logger (so it's not silently missing)
  missing_sites <- setdiff(sort(unique(scored$sitekey)), best$sitekey)
  if (length(missing_sites) > 0)
    message("  [ref-finder] NO qualifying reference logger for: ",
            paste(missing_sites, collapse = ", "),
            " (relax depth_sd_max / min_hours, or accept satellite Kd there)")
  
  cat("\n── RECOMMENDED reference logger per site (CONSTANT — from output) ──\n")
  print(as.data.frame(best %>%
                        transmute(sitekey, hoboID_env,
                                  shadow     = round(shadow_mean, 3),
                                  shadow_min = round(shadow_min,  3),
                                  ambient    = round(ambient_mean, 3),
                                  depth_m    = round(depth_mean,  2),
                                  n_hours,
                                  score      = round(open_score,  3))),
        row.names = FALSE)
  
  # ready-to-paste named vector for REFERENCE_LOGGER
  cat("\n# Paste into REFERENCE_LOGGER (stable across re-runs / toggle flips):\n")
  cat("REFERENCE_LOGGER <- c(\n")
  if (nrow(best) > 0) {
    cat(paste0("  ", best$sitekey, ' = "', best$hoboID_env, '"',
               c(rep(",", nrow(best) - 1), ""), collapse = "\n"))
  }
  cat("\n)\n")
  
  invisible(best)
}

# Score from `output` (the full, undropped prediction table) so the
# recommendation is identical whether or not USE_REFERENCE_KD is on.
ref_candidates <- find_reference_loggers(output)



# ══════════════════════════════════════════════════════════════════════════════
#  PER-SITE FIGURE 1 — measured vs predicted PAR per logger
#  One figure PER site/year, saved to file. Loops over dem_store keys.
# ══════════════════════════════════════════════════════════════════════════════
plot_par_panels_one_site <- function(key, val_all) {
  val <- val_all %>% filter(sitekey == key)
  if (nrow(val) == 0) { message("  (no validation rows for ", key, ")"); return(invisible()) }
  
  date_range  <- format(range(val$datetime_pos), "%Y-%m-%d %H:%M")
  out_px      <- max(dim(dem_store[[key]]$dem$matrix))
  kd_disp <- if (USE_REFERENCE_KD) "per-hour"
  else ifelse(key %in% names(KDPAR_BY_KEY), round(KDPAR_BY_KEY[[key]],3), KD_CONSTANT)
  plot_title  <- paste0(
    key, " | ", kd_method_label(), " | ", date_range[1], " to ", date_range[2], "\n",
    "Kd=", kd_disp, " | zscale=", round(dem_store[[key]]$dem$zscale, 4), " (UTM)",
    " | interval=", TIMESTEP_INTERVAL, " | ~", out_px, "px (target ", TARGET_PIXELS, ")",
    " | maxsearch=", MAXSEARCH_FACT, " | w_dir=", W_DIRECT, " w_amb=", W_AMBIENT
  )
  
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
  
  png(paste0("par_panels_", key, ".png"),
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


plot_par_panels_one_site <- function(key, val_all) {
  val <- val_all %>% filter(sitekey == key)
  if (nrow(val) == 0) { message("  (no validation rows for ", key, ")"); return(invisible()) }
  
  date_range  <- format(range(val$datetime_pos), "%Y-%m-%d %H:%M")
  out_px      <- max(dim(dem_store[[key]]$dem$matrix))
  kd_disp <- if (USE_REFERENCE_KD) "per-hour"
  else ifelse(key %in% names(KDPAR_BY_KEY), round(KDPAR_BY_KEY[[key]],3), KD_CONSTANT)
  plot_title  <- paste0(
    key, " | ", kd_method_label(), " | ", date_range[1], " to ", date_range[2], "  |  ",
    "Kd=", kd_disp,
    " | zscale=", round(dem_store[[key]]$dem$zscale, 4),
    " | maxsearch=", MAXSEARCH_FACT, " | w_dir=", W_DIRECT, " w_amb=", W_AMBIENT
  )
  
  ylim_shared <- range(c(val$pred_par, val$light_mol), na.rm = TRUE)
  hobos       <- sort(unique(val$hoboID_env))
  n_hobos     <- length(hobos)
  
  # ── landscape grid: wide, few rows ──────────────────────────────────────────
  n_cols      <- min(n_hobos, 8)                 # up to 8 panels across
  n_rows      <- ceiling(n_hobos / n_cols)
  legend_cell <- n_hobos + 2
  plot_cells  <- c(2:(n_hobos + 1), rep(0, n_cols * n_rows - n_hobos))
  
  layout_mat <- rbind(
    rep(1, n_cols),                              # title strip
    matrix(plot_cells, nrow = n_rows, ncol = n_cols, byrow = TRUE),
    rep(legend_cell, n_cols)                     # legend strip
  )
  
  # ── 16:9 canvas (1600x900 at res 150 ≈ a full slide) ────────────────────────
  png(paste0("par_panels_", key, ".png"),
      width = 1600, height = 900, res = 150)
  layout(layout_mat, heights = c(0.5, rep(3, n_rows), 0.5))
  
  par(mar = c(0, 0, 0, 0)); plot.new()
  text(0.5, 0.5, plot_title, cex = 1.0, font = 2, adj = 0.5)
  
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


# ══════════════════════════════════════════════════════════════════════════════
#  PER-SITE FIGURE 2 — shadow maps + per-timestep bar charts
#  One figure PER site/year, saved to file. Reuses cached DEM/hobo from dem_store.
# ══════════════════════════════════════════════════════════════════════════════
plot_shadow_multipanel_one_site <- function(key, output_all, val_all) {
  e          <- dem_store[[key]]
  dem        <- e$dem
  hobo_sf    <- e$hobo_sf
  hobo_wgs84 <- e$hobo_wgs84
  
  out <- output_all %>% filter(sitekey == key)
  val <- val_all    %>% filter(sitekey == key)
  if (nrow(out) == 0) { message("  (no output rows for ", key, ")"); return(invisible()) }
  
  sun_steps <- out %>% distinct(datetime_pos, sun_altitude, sun_azimuth) %>%
    arrange(datetime_pos)
  n_steps <- nrow(sun_steps)
  n_cols  <- 4
  n_rows  <- ceiling(n_steps / n_cols)
  
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
         paste0(key, "  |  ", kd_method_label(),
                "  |  Kd=", if (USE_REFERENCE_KD) "per-hour"
                else ifelse(key %in% names(KDPAR_BY_KEY), round(KDPAR_BY_KEY[[key]],3), KD_CONSTANT),
                "  | zscale=", round(dem$zscale, 4), " (UTM)",
                "  |  ~", max(dim(dem$matrix)), "px (target ", TARGET_PIXELS, ")",
                "  |  maxsearch=", MAXSEARCH_FACT,
                "  |  w_dir=", W_DIRECT, " w_amb=", W_AMBIENT),
         cex = 1.2, font = 2, adj = 0.5)
    
    for (i in seq_len(n_steps)) {
      ts      <- sun_steps[i, ]
      alt_deg <- ts$sun_altitude * 180 / pi
      az_deg  <- ts$sun_azimuth  * 180 / pi
      
      par(mar = c(0, 1, 5, 1))
      shadow <- ray_shade(dem$matrix, sunaltitude = alt_deg, sunangle = az_deg,
                          zscale = dem$zscale, maxsearch = MAXSEARCH_FACT)
      shadow[!dem$reef_mask] <- NA
      shadow_rast <- rast(shadow)
      ext(shadow_rast) <- ext(dem$rast)
      crs(shadow_rast) <- crs(dem$rast)
      
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
  
  png(paste0("shadow_multipanel_", key, ".png"),
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
#  SHAPE = site/year, COLOUR = depth, every sensor labelled.
# ══════════════════════════════════════════════════════════════════════════════
# ── COLOUR TOGGLE: "depth" or "shadow" ────────────────────────────────────────
COLOR_BY <- "depth"   # "depth" = water column (m); "shadow" = mean shadow fraction

# ── DLI UNIT / INTEGRATION-WINDOW FIX ─────────────────────────────────────────
#  OLD BUG: measured DLI summed the raw 5-min light_mol from env_long_test
#  (x 300 s) while predicted DLI summed hourly pred_par (x 3600 s). The two
#  therefore integrated over DIFFERENT row sets and DIFFERENT time bases:
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
#  Only hours present in `validation` (sun_alt>15 AND a finite prediction) count
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

# SHAPE per site/year
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

dev.new(width = 8, height = 7)
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
                   " | maxsearch=", MAXSEARCH_FACT,
                   " | w_dir=", W_DIRECT, " w_amb=", W_AMBIENT),
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

# ══════════════════════════════════════════════════════════════════════════════
#  POOLED FIGURE 2 — observed vs predicted PAR at selected hours
#  SHAPE = site/year, COLOUR = depth, every sensor labelled.
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
  
  dev.new(width = 8, height = 7)
  par(mar = c(5, 5, 4, 7))
  plot(time_df$light_mol, time_df$pred_par,
       pch = pch_map2[time_df$sitekey], col = col2, cex = 1.4, lwd = 2,
       xlim = lim2, ylim = lim2,
       xlab = "Observed PAR (µmol)", ylab = "Predicted PAR (µmol)",
       main = paste0("PAR at 8/10/12/15/17h — ", length(keys2),
                     " site/years (colour = hour) | ", kd_method_label()),
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
  # hour legend (colour) — the new primary key
  legend(lim2[2] + diff(lim2) * 0.04, lim2[2],
         legend = paste0(hour_levels, ":00"),
         col = hour_cols[as.character(hour_levels)],
         pch = 19, pt.cex = 1.4, bty = "n", title = "Hour")
  # site legend (shape)
  legend(lim2[2] + diff(lim2) * 0.04, lim2[2] - diff(lim2) * 0.28,
         legend = keys2, pch = pch_map2[keys2], pt.cex = 1.3, pt.lwd = 2,
         bty = "n", title = "Site/year")
  par(xpd = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
#  OPTIONAL — FIT W_DIRECT / W_AMBIENT (and Kd). Pools ALL site/years in
#  `validation`, so multi-site just gives it a bigger fitting set. Unchanged.
# ══════════════════════════════════════════════════════════════════════════════
fit_ambient_weights <- function(validation,
                                kd_grid = seq(0.05, 0.60, by = 0.01),
                                wamb_grid = seq(0, 0.95, by = 0.02)) {
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

# usage after pipeline:
#   fit <- fit_ambient_weights(validation)
#   W_DIRECT  <- 1 - fit$best_wamb
#   W_AMBIENT <- fit$best_wamb
#   # then re-run run_pipeline() to apply the fitted weights