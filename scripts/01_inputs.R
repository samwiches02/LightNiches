# ══════════════════════════════════════════════════════════════════════════════
#  01_inputs.R — build the hourly input table
#
#  Reads : data/env_long_depth_kd_longlat.csv, data/tides/, data/kdpar_lookup.csv,
#          data/kd_raw_per_wavelength.csv, data/geom_new.csv
#  Writes: output/hourly_inputs.rds   one row per logger per hour
#          output/capture_tide.csv    tide at each photogrammetry capture
#          output/kd_satellite.csv    satellite KdPAR per sampling event
#          output/kd_bands.rds        band-by-band Kd (supplement)
#          output/geom.rds            logger positions
#
#  [NEW] as a script. The steps are yours (they ran at the top of
#  master_script.R and inside run_pipeline); the functions they call are tagged
#  in R/functions_inputs.R.
#
#  PROVENANCE TAGS:  [YOURS] your code from master_script.R, logic unchanged
#                    [YOURS, EDITED] your code with listed changes   [NEW] Claude
#  Comments with no "[Claude]" marker are yours, copied from master_script.R.
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_inputs.R"))
suppressPackageStartupMessages({ library(dplyr); library(lubridate) })

stage_banner("01  INPUTS")

# ── Tides and the capture-tide correction ────────────────────────────────────
# [YOURS] — your "Fitting harmonic tide models" and capture-tide lookup steps
# LOOKUP TABLE — one capture_tide_m per site/year. Printed so you can eyeball it before runs.
tide_fit     <- fit_tide_models_from_csv()
capture_tide <- build_capture_tide(tide_fit)
print(capture_tide[capture_tide$sitekey %in% SITES, ], row.names = FALSE)
write.csv(capture_tide, out_path("capture_tide.csv"), row.names = FALSE)

# ── Kd ───────────────────────────────────────────────────────────────────────
# [YOURS, EDITED] — your per-site KdPAR printout, now saved to a file
kd_sat <- read_kd_satellite()
write.csv(data.frame(sitekey = names(kd_sat), kd_sat = kd_sat, column = KDPAR_COL),
          out_path("kd_satellite.csv"), row.names = FALSE)
saveRDS(read_kd_bands(), out_path("kd_bands.rds"))

# ── Hourly logger table ──────────────────────────────────────────────────────
# [YOURS, EDITED] — your hourly schedule plus the light_mol averaging from your
#   validation join, done in one place
env    <- read_env_long()
hourly <- build_hourly_inputs(env, capture_tide, kd_sat)

message("\nHourly rows per sampling event:")
print(table(hourly$sitekey))
message("Rows with water depth <= 0 (logger exposed): ", sum(hourly$water_depth_m <= 0, na.rm = TRUE))

saveRDS(hourly, out_path("hourly_inputs.rds"))

# ── Logger positions ─────────────────────────────────────────────────────────
# [YOURS, EDITED] geom loading. [NEW] the check for loggers with no position
geom <- read_geom()
missing_pos <- hourly %>% distinct(sitekey, hoboID_env) %>%
  anti_join(geom, by = c("sitekey", "hoboID_env" = "hoboID"))
if (nrow(missing_pos) > 0) {
  warning(nrow(missing_pos), " logger(s) have light data but no position in geom_new.csv; ",
          "they will be dropped in 02. See output/loggers_without_position.csv")
  write.csv(missing_pos, out_path("loggers_without_position.csv"), row.names = FALSE)
}
saveRDS(geom, out_path("geom.rds"))

write_manifest("01_inputs", extra = c(
  sprintf("hourly rows: %d across %d events", nrow(hourly), length(unique(hourly$sitekey))),
  sprintf("loggers without positions: %d", nrow(missing_pos))))
message("01 done.")
