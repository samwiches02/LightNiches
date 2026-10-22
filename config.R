# ══════════════════════════════════════════════════════════════════════════════
#  config.R — EVERY setting for the run lives here, and only here.
#
#  PROVENANCE TAGS (every section has one):
#    [YOURS]          your values and comments from master_script.R
#    [YOURS, EDITED]  your values, with the listed changes
#    [NEW]            added by Claude
#  Comments with no "[Claude]" marker are yours, copied from master_script.R.
#
#  [Claude] Rule: no script or function assigns any of these names. In
#  master_script.R the weighting constants were assigned twice, ~300 lines
#  apart, so the values in the "RUN" block were not the ones that ran.
#
#  [Claude] Removed: USE_CONSTANT_KD / USE_REFERENCE_KD / USE_SPECTRAL_KD /
#  USE_AMBIENT toggles (every model now runs in one pass), KDPAR_OVERRIDE,
#  USE_ALT_DEPENDENT_WEIGHTS, USE_DEPTH_DEPENDENT_WEIGHTS, W_AMBIENT_MIN/MAX/P,
#  TAU_SCALE, ALT_EXP, DEPTH_EXP, W_DIRECT/W_AMBIENT, N_REF_LOGGERS,
#  REFERENCE_LOGGER (hand-pasted list) and REF_AGG_FUN (now always the median).
# ══════════════════════════════════════════════════════════════════════════════


# ── Sampling events ───────────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] your `wanted` table, as keys
# [Claude, 20 Oct 2026] every event with light data, a DEM and a capture time
SITES <- c("CB18", "CB19",
           "GT17", "GT19",
           "HS17", "HS19",
           "L117", "L119",
           "N317", "N318", "N319",
           "RS18",
           "SE17", "SE18",
           "TB18", "TB19",
           "TM19")

# [NEW] Every event with data that is NOT in SITES needs a written reason here.
#   03_models.R prints this table.
EXCLUDED_EVENTS <- c(
  GT18 = "no DEM",
  HS18 = "no DEM",
  L118 = "no DEM",
  CB17 = "no photogrammetry capture time (needed for the capture-tide depth correction)",
  RS17 = "no photogrammetry capture time (needed for the capture-tide depth correction)",
  TB17 = "no photogrammetry capture time (needed for the capture-tide depth correction)",
  TM17 = "no photogrammetry capture time (needed for the capture-tide depth correction)",
  SE19 = "no photogrammetry capture time (needed for the capture-tide depth correction)",
  RS19 = "no photogrammetry capture time; also only 6 light records per logger",
  TM18 = paste("no reference logger: every logger passing the geometry rules is",
               "shallower than REF_DEPTH_MIN_M (1 m) at low tide (0.6-0.8 m);",
               "decided before any metric was computed")
)

# ── Paths ─────────────────────────────────────────────────────────────────────
# [NEW] relative to the project root (here::here() finds it via the .here file).
#   Replaces setwd("~/Desktop/...") and the /Users/samueleszenyi/... paths.
DATA_DIR  <- here::here("data")
OUT_DIR   <- here::here("output")

ENV_LONG_CSV     <- file.path(DATA_DIR, "env_long_depth_kd_longlat.csv")
GEOM_CSV         <- file.path(DATA_DIR, "geom_new.csv")
KDPAR_LOOKUP_CSV <- file.path(DATA_DIR, "kdpar_lookup.csv")
KD_BANDS_CSV     <- file.path(DATA_DIR, "kd_raw_per_wavelength.csv")
# Directory holding tides_Lizard{2017,2018,2019}.csv
TIDE_CSV_DIR     <- file.path(DATA_DIR, "tides")
#create filepath for DEMS   [Claude: files named <SITE><yy>.tif, e.g. CB19.tif]
DEM_DIR          <- file.path(DATA_DIR, "dems")

LOCAL_TZ <- "Australia/Brisbane"


# ── Light sensors ─────────────────────────────────────────────────────────────
# [NEW] — the value is from read_environment.R:
#   light$mol_light <- 0.0185*light$light # translation coefficient from Tilman
#   [Claude] The thesis cites Howell et al. (1983) for this. Make the code,
#   this comment and the paper cite the SAME source.
LUX_TO_PAR <- 0.0185


# ── Time and sun ──────────────────────────────────────────────────────────────
# [YOURS, EDITED]
#   - [Claude] the 20° cutoff was hard-coded in run_pipeline(); now evaluated
#     at mid-hour
TIMESTEP_INTERVAL <- "1 hour"    # the timestep interval to measure
# discard hours where the sun is less than 20 degrees over the horizon.
# at super low angles, the shadows hit DEM imperfections harder (long
# shadows from tiny noise)
SUN_ALT_MIN_DEG   <- 20


# ── Optics ────────────────────────────────────────────────────────────────────
# [NEW]
N_WATER     <- 1.34    # [Claude] refractive index of seawater (Snell's law)
USE_FRESNEL <- FALSE   # [Claude] TRUE = multiply by air-water transmittance; applied
                       # to predictions AND the reference-Kd inversion together

# [YOURS]
KD_CONSTANT <- 0.2     # value used for the constant-Kd model (Macdonald 2015)
# which KdPAR variant to apply:
#   "kdpar_photon" = photon-weighted full-band integral
#   "kdpar_mean"   = unweighted mean of the 10 bands
#   "kdpar_geom"   = geometric mean, red-outlier-resistant
#   "kdpar_490corr" = 490 correlated band
KDPAR_COL   <- "kdpar_photon"


# ── DEM / ray-tracing geometry ────────────────────────────────────────────────
# [YOURS] — values unchanged
# all DEMs are reprojected to UTM zone 55S before processing
# EPSG:32755 = WGS84 UTM Zone 55S, correct for Lizard Island (~145°E, ~15°S)
UTM_CRS               <- "EPSG:32755"

# [NEW] Claude. CRS to ASSIGN to a DEM before reprojection, when the file's own
#   CRS is missing or wrong. Keys are an event ("CB17") or a year ("2017");
#   an event entry wins over its year. Leave empty if the files carry the right
#   CRS; they are reprojected to UTM_CRS automatically either way. Example:
#   DEM_CRS_OVERRIDE <- list("2017" = "EPSG:28355")   # GDA94 / MGA zone 55
# 2017 batch: checked 20 Oct 2026 from the GeoTIFF keys of all eight 2017 files.
#   They are in the local WGS84 Transverse Mercator (origin 14.69897 S,
#   145.44826 E, k = 1, no false easting/northing): the same frame as the logger
#   positions in geom_new.csv (TMERC_CRS). Their CRS is user-defined in the
#   file, so it is set explicitly here; the value is identical to the file's.
#   2018/2019 files are WGS84 lon/lat (EPSG:4326) and need nothing.
DEM_CRS_OVERRIDE <- list("2017" = paste0("+proj=tmerc +lat_0=-14.6989736557007",
                                         " +lon_0=145.448257446289 +k=1 +x_0=0 +y_0=0",
                                         " +datum=WGS84 +units=m +no_defs"))
# Sanity checks after reprojection (standardize_dem_crs)
LIZARD_BBOX      <- c(lon_min = 145.35, lon_max = 145.55, lat_min = -14.80, lat_max = -14.55)
DEM_MAX_RELIEF_M <- 20

#  REEF_RADIUS_M     : every site gets a 6.5 m radius (~13 m diameter) disc
#                      centred on the data centroid. Identical physical reef
#                      area on every site, regardless of how the DEM was
#                      processed.
REEF_RADIUS_M         <- 6.5
CROP_MARGIN_M         <- 0.5   # cosmetic buffer around the 6.5 m disc

# DOWNSAMPLING: target the REEF DISC's pixel resolution, not the matrix.
# Set REEF_RADIUS_PX_TARGET (px on the downsampled grid that the 6.5 m reef
# should fill), and derive the downsample factor per-DEM from the NATIVE pixel
# size so the disc ends up at that resolution on every site.
# Bigger value = finer disc = slower raytrace (cost scales ~N^2).
REEF_RADIUS_PX_TARGET <- 800   # target pixel radius of the reef disc
MIN_DOWNSAMPLE        <- 1     # never upsample; 1 = keep native if already finer

#  MAXSEARCH_M       : direct-beam rays travel the same physical distance
#                      per site. Converted to pixels per-DEM via zscale
#                      inside make_shadow().
MAXSEARCH_M           <- 13
#  AMBIENT_MAXSEARCH_M: same for the ambient layer; converted in make_ambient().
AMBIENT_MAXSEARCH_M   <- 6.5
AMBIENT_SUNBREAKS     <- 30    # rays sent in a circle around each pixel (rayshader default is 24)

#Buffer to extract shadow and ambient values (always matching) in meters
BUFFER_M         <- 0.06
# [NEW] extra buffers stored for the sensitivity check
BUFFER_SENS_M    <- c(0.02, 0.10)

# [NEW] Also ray-trace with the UNrefracted (in-air) sun so the paper can show
#   what refraction changes (model M3r_unrefracted). Doubles ray-tracing time.
RUN_UNREFRACTED_DIAGNOSTIC <- TRUE

# [NEW] Re-ray-trace a site even if it is already saved in output/raytrace/
FORCE_RAYTRACE <- FALSE

# [YOURS, EDITED] — your LOGGER_OFFSETS tribble, as a data.frame
#   apply to lon/lat (geom is in local tmerc, units = metres, so dx/dy add directly)
#   [Claude] document every row you add in the Methods.
LOGGER_OFFSETS <- data.frame(site = character(), year = integer(),
                             hoboID = character(), dx_m = numeric(), dy_m = numeric())


# ── Tides and the DEM capture-tide correction ────────────────────────────────
# [YOURS]
# [Claude, 20 Oct 2026] was "hc37". Each tide record is only 16-26 days long,
#   too short to separate 37 constituents (Rayleigh criterion): the hc37 fit
#   matched the record to <2 cm but gave 4-6 m and negative tides 2-4 days
#   outside it (all 2017 captures; HS19/CB19/N319/L119/TB19 are 2-4 days out
#   too). Holding out the first 4 days of each record and predicting them:
#   these 8 constituents give RMSE 15 cm (2017), 11 cm (2018), 8 cm (2019);
#   the 37-constituent fit gives errors of metres. Can be a TideHarmonics set
#   name ("hc7", "hc37") or a vector of constituent names.
TIDE_CONSTITS    <- c("M2", "S2", "N2", "K1", "O1", "Q1", "M4", "MS4")
EXTRAP_WARN_DAYS <- 7        # warn if a capture is >this far outside the fit window

# [YOURS, EDITED] — your photogrammetry tribble
#   - [Claude] keyed by sitekey only; the 2017 events and L119 (not in the
#     current analysis) were dropped; add them back if needed
# ── PHOTOGRAMMETRY CAPTURE TIMES (local Brisbane time) ────────────────────────
PHOTOGRAMMETRY <- data.frame(
  sitekey    = c("TB18", "CB18", "RS18", "TB19", "CB19", "GT19", "HS19",
                 "N318", "N319", "TM18", "TM19",
                 # [Claude] restored from master_script.R's tribble
                 "GT17", "HS17", "L117", "L119", "N317", "SE17", "SE18"),
  capture_dt = c("2018-11-30 14:39", "2018-11-02 08:32", "2018-12-02 16:15",
                 "2019-11-03 13:54", "2019-11-03 09:09", "2019-11-08 10:07",
                 "2019-11-02 10:24", "2018-11-28 13:09", "2019-11-03 10:24",
                 "2018-11-23 09:24", "2019-11-12 10:43",
                 "2017-11-08 11:22", "2017-11-08 16:43", "2017-11-08 16:08",
                 "2019-11-03 16:20", "2017-11-07 10:01", "2017-11-07 09:13",
                 "2018-11-25 08:55"),
  stringsAsFactors = FALSE
)

# [YOURS]
#  Keyed by `abbrev`. Any abbrev NOT listed keeps its harmonic value.
CAPTURE_TIDE_OVERRIDE <- c(
  CB18 = 1.13   # <-- REQUIRED for CB18: external Lizard tide @ 2018-11-02 08:32
                # [Claude] TODO: record the source of this value
)


# ── Reference loggers ─────────────────────────────────────────────────────────
# [NEW] — replaces the hand-pasted REFERENCE_LOGGER list and the
#   find_reference_loggers() weights. Geometry only; set these BEFORE looking at
#   any model metrics. 03_models.R builds the pool, saves it to
#   output/reference_pool.csv and STOPS the first time, so you can check it.
#   Once the file exists it is read, never rebuilt.
REF_POOL_MAX    <- 8      # candidates kept per site (used for the logger-number test)
N_REF           <- 5      # loggers used for the main reference-Kd model (your N_REF_LOGGERS)
REF_MIN_HOURS   <- 5      # your min_hours gate
REF_MIN_MIDDAY  <- 2      # your >= 2 midday hours gate (hours 10:00-14:00)
REF_AMBIENT_MIN <- 0.90   # mean ambient (sky-view) fraction in the buffer
REF_SHADOW_MIN  <- 0.80   # refracted direct-light fraction in the logger's WORST hour
REF_DEPTH_MIN_M <- 1.0    # shallowest water depth during deployment


# ── Model fitting and statistics ──────────────────────────────────────────────
# [NEW]
W_GRID      <- seq(0, 1, by = 0.05)   # ambient weight, chosen by leave-one-event-out CV
BOOT_B      <- 1000                   # cluster-bootstrap replicates (loggers within events)
REF_KS      <- c(1, 2, 3, 5)          # reference-logger numbers to test
REF_N_DRAW  <- 200                    # random draws per k
SEED        <- 2026
RUN_SPECTRAL_SUPPLEMENT <- TRUE       # band-by-band Kd, reported in the supplement only


# ── Figures ───────────────────────────────────────────────────────────────────
# [YOURS, EDITED] — your plot_shadow_3hr_with_zoom("CB19", ...) choices
FIG_SITE        <- "CB19"         # event shown in the shadow-map figure
FIG_HOURS       <- c(8, 12, 16)   # local hours shown in the shadow-map figure
# [NEW]
FIG_MAIN_MODEL  <- "M4r"          # model used for residual and anomaly figures



# ── Test run switch ───────────────────────────────────────────────────────────
# TRUE = quick CB-only run into output_test/; FALSE = full run into output/
TEST_RUN <- FALSE
if (TEST_RUN) {
  SITES   <- c("CB18", "CB19")
  EXCLUDED_EVENTS <- c(EXCLUDED_EVENTS, TEST = "test run: CB only")
  OUT_DIR <- here::here("output_test")   # keeps test files away from output/
  REEF_RADIUS_PX_TARGET <- 300           # coarser DEM: ray-tracing ~7x faster
  BOOT_B     <- 100                      # fewer bootstrap replicates
  REF_N_DRAW <- 20                       # fewer logger-number draws
  REF_KS <- c(1, 2, 3)
}


# ── Run labels and sensitivity runs ───────────────────────────────────────────
# [NEW] Claude, 22 Oct 2026
#  RUN_LABEL "main" = the pre-registered analysis, written to output/.
#  Any other label writes stages 03, 04 and 07 to output/runs/<label>/ and
#  READS the stage 01/02 products (cache, DEMs, reference pool, eval rows) from
#  output/, so nothing in output/ is overwritten. Stages 01/02 refuse to run
#  under a label.
#  W_AMBIENT_FIXED: NA = ambient weight chosen by leave-one-event-out CV (as
#  pre-registered). A number fixes it for every event (exploratory). Setting it
#  without a label gives the label "w_ambient_<value>" automatically.
#  Both can be set for one run from the command line, without editing this file:
#     Rscript run_all.R 03 07 W_AMBIENT_FIXED=0.7
RUN_LABEL       <- "main"
W_AMBIENT_FIXED <- NA

for (.nm in c("RUN_LABEL", "W_AMBIENT_FIXED")) {          # command-line overrides
  .v <- Sys.getenv(paste0("LN_", .nm), "")
  if (nzchar(.v)) assign(.nm, utils::type.convert(.v, as.is = TRUE))
}
rm(list = intersect(c(".nm", ".v"), ls(all.names = TRUE)))
if (is.finite(W_AMBIENT_FIXED)) {
  if (RUN_LABEL == "main") RUN_LABEL <- paste0("w_ambient_", W_AMBIENT_FIXED)
  W_GRID <- W_AMBIENT_FIXED
}
BASE_OUT_DIR <- OUT_DIR                                     # stage 01/02 products live here
if (RUN_LABEL != "main") OUT_DIR <- file.path(BASE_OUT_DIR, "runs", RUN_LABEL)
RUN_ID <- Sys.getenv("LN_RUN_ID", format(Sys.time(), "%Y%m%d_%H%M%S"))  # shared by one run_all.R call
