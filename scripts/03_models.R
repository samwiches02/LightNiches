# ══════════════════════════════════════════════════════════════════════════════
#  03_models.R — reference loggers, model ladder, cross-validation, metrics
#
#  Reads : output/cache.rds, output/kd_bands.rds
#  Writes: output/reference_pool.csv   FROZEN geometry-only pool (built once)
#          output/eval_keys.csv        FROZEN evaluation rows (built once)
#          output/predictions.rds      every model's prediction on the same rows
#          output/metrics_*.csv        pooled, per-event, daytime-integrated,
#                                      within-site, bootstrap CIs, events improved
#          output/ref_number.csv       error vs number of reference loggers
#
#  The first run STOPS after writing reference_pool.csv, before any metric
#  exists. Check the pool, then run again. To rebuild it you must delete the
#  file yourself, which keeps the choice deliberate and visible in git.
#
#  [NEW] overall: written by Claude. It replaces your toggle-and-rerun workflow
#  (one run per Kd method / weighting, logged to metrics_log.csv) with one pass
#  that scores every model on the same rows. Pieces that are yours are tagged.
#
#  PROVENANCE TAGS:  [YOURS] your code from master_script.R, logic unchanged
#                    [YOURS, EDITED] your code with listed changes   [NEW] Claude
#  Comments with no "[Claude]" marker are yours, copied from master_script.R.
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "functions_models.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

stage_banner("03  MODELS")
message("Run label: ", RUN_LABEL, "   output folder: ", OUT_DIR,
        "\nAmbient weight: ", if (is.finite(W_AMBIENT_FIXED))
          paste("FIXED at", W_AMBIENT_FIXED, "(exploratory)") else "leave-one-event-out CV")

cache <- readRDS(in_path("cache.rds"))
cache$T_surf <- if (USE_FRESNEL) fresnel_T(cache$sun_altitude) else 1

message("Excluded sampling events and reasons:")
print(EXCLUDED_EVENTS)

# ── 1. Reference pool (geometry only, frozen) ────────────────────────────────
# [NEW] replaces your BOOTSTRAP WORKFLOW (run with satellite Kd, paste the
#   finder's output into REFERENCE_LOGGER, re-run)
pool_file <- in_path("reference_pool.csv")
if (!file.exists(pool_file)) {
  pool <- select_reference_pool(cache)
  write.csv(pool, pool_file, row.names = FALSE)
  message("\nReference pool written to ", pool_file, ":")
  print(as.data.frame(pool %>% count(sitekey, name = "n_candidates")))
  none <- setdiff(SITES, pool$sitekey)
  if (length(none) > 0)
    message("No candidate passed the thresholds at: ", paste(none, collapse = ", "),
            ". If you relax a threshold in config.R, do it now, delete the file and",
            " re-run, and report the final thresholds in the Methods.")
  end_stage(paste0("\nSTOPPING before any metric is computed. Check the pool (and the",
                   " shadow maps from 04_figures.R if useful), then run 03 again."),
            STATUS_PAUSE)
}
pool <- read.csv(pool_file, stringsAsFactors = FALSE)
refs <- pool %>% filter(rank <= N_REF) %>% select(sitekey, hoboID_env)
message("Reference loggers used for Kd: ", nrow(refs), " across ",
        length(unique(refs$sitekey)), " events")

# ── 2. Reference Kd (hourly, and a site-constant version) ────────────────────
# [YOURS, EDITED] — your reference-Kd method; see reference_kd()
kd <- reference_kd(cache, refs)
message("Hours filled with the site-median reference Kd: ", sum(kd$kd_filled), " of ", nrow(kd))
cache <- cache %>%
  left_join(kd %>% select(sitekey, datetime_pos, kd_ref, kd_ref_site, kd_filled),
            by = c("sitekey", "datetime_pos"))

# ── 3. Evaluation rows: identical for every model (frozen) ───────────────────
# [NEW] — your code dropped reference loggers only in reference-Kd runs, so
#   satellite and reference runs were scored on different rows
#  Every pool logger is excluded (not only the N_REF used), so the logger-number
#  experiment and the main models are scored on loggers never eligible as references.
sens_sfx <- buffer_suffix(BUFFER_SENS_M)
needed <- c("surface_par", "light_mol", "water_depth_m", "shadow_refr", "ambient",
            "kd_ref", "kd_ref_site",          # [Claude] kd_sat not required (NA for 2017)
            paste0("shadow_refr", sens_sfx), paste0("ambient", sens_sfx))
if (RUN_UNREFRACTED_DIAGNOSTIC) needed <- c(needed, "shadow_air")
eval_d <- cache %>%
  anti_join(pool, by = c("sitekey", "hoboID_env")) %>%
  filter(if_all(all_of(needed), is.finite), water_depth_m > 0)

keys_file <- in_path("eval_keys.csv")
eval_keys <- eval_d %>% distinct(sitekey, hoboID_env, datetime_pos) %>%
  arrange(sitekey, hoboID_env, datetime_pos)
if (!file.exists(keys_file)) {
  write.csv(transform(eval_keys, datetime_pos = format(datetime_pos, "%Y-%m-%d %H:%M:%S")),
            keys_file, row.names = FALSE)
} else {
  frozen <- read.csv(keys_file, stringsAsFactors = FALSE)
  if (nrow(frozen) != nrow(eval_keys))
    stop("Evaluation rows changed (", nrow(frozen), " frozen vs ", nrow(eval_keys),
         " now). Something upstream changed; find out what before deleting eval_keys.csv.")
}
N_EVAL <- nrow(eval_d)
message("Evaluation set: ", N_EVAL, " logger-hours from ",
        nrow(distinct(eval_d, sitekey, hoboID_env)), " loggers in ",
        length(unique(eval_d$sitekey)), " events")

# ── 4. Model ladder: every model predicts the same rows ──────────────────────
# [NEW] — the equation is yours (predict_par); the ladder, the single
#   cross-validated weight and the diagnostics are new
add <- function(d, model, pred, w_fold = NA_real_) {
  stopifnot(length(pred) == nrow(d), all(is.finite(pred)))
  d %>% transmute(model = model, sitekey, hoboID_env, datetime_pos,
                  light_mol, water_depth_m, pred = pred, w_fold = w_fold)
}
loeo_add <- function(model, kd_col, sfx = "") {
  r <- loeo_w(eval_d, kd_col, paste0("shadow_refr", sfx), paste0("ambient", sfx))
  r <- r[match(paste(eval_d$sitekey, eval_d$hoboID_env, eval_d$datetime_pos),
               paste(r$sitekey, r$hoboID_env, r$datetime_pos)), ]
  add(eval_d, model, r$pred, r$w_fold)
}

# [Claude] satellite models: only events with a satellite Kd (none for 2017)
eval_s <- eval_d %>% filter(is.finite(kd_sat))
message("Satellite models scored on ", nrow(eval_s), " of ", N_EVAL, " logger-hours (",
        length(unique(eval_s$sitekey)), " of ", length(unique(eval_d$sitekey)),
        " events have a satellite Kd)")
r_s <- loeo_w(eval_s, "kd_sat")
r_s <- r_s[match(paste(eval_s$sitekey, eval_s$hoboID_env, eval_s$datetime_pos),
                 paste(r_s$sitekey, r_s$hoboID_env, r_s$datetime_pos)), ]

preds <- list(
  add(eval_d, "M0",  predict_par(eval_d, NULL,       "none")),
  add(eval_d, "M1",  predict_par(eval_d, "kd_const", "none")),
  # [NEW] 29 Oct 2026: BASELINE model = constant Kd + direct shading only (like the
  #   thesis baseline). M1d -> M3r isolates the Kd source; M1 -> M1d the shading.
  add(eval_d, "M1d", predict_par(eval_d, "kd_const", "direct")),
  add(eval_s, "M2s", predict_par(eval_s, "kd_sat",   "none")),
  add(eval_d, "M2r", predict_par(eval_d, "kd_ref",   "none")),
  add(eval_s, "M3s", predict_par(eval_s, "kd_sat",   "direct")),
  add(eval_d, "M3r", predict_par(eval_d, "kd_ref",   "direct")),
  add(eval_s, "M4s", r_s$pred, r_s$w_fold),
  loeo_add("M4r", "kd_ref"),
  loeo_add("M4r_siteKd", "kd_ref_site")                    # hourly vs site-constant Kd
)
if (RUN_UNREFRACTED_DIAGNOSTIC)
  preds <- c(preds, list(add(eval_d, "M3r_unrefracted",
                             predict_par(eval_d, "kd_ref", "direct_air"))))
for (b in BUFFER_SENS_M)                                    # buffer sensitivity
  preds <- c(preds, list(loeo_add(paste0("M4r_buffer", buffer_suffix(b)), "kd_ref",
                                  buffer_suffix(b))))
if (RUN_SPECTRAL_SUPPLEMENT) {
  kd_bands <- readRDS(in_path("kd_bands.rds"))
  sp <- eval_d$surface_par * eval_d$T_surf *
    spectral_atten(eval_d$water_depth_m, eval_d$sitekey, kd_bands)
  if (all(is.finite(sp))) preds <- c(preds, list(add(eval_d, "S_spectralKd", sp)))
  else message("Spectral supplement skipped: some events lack band Kd.")
}
preds <- bind_rows(preds)
saveRDS(preds, out_path("predictions.rds"))

# [NEW] Claude, 22 Oct 2026 — dated validation sheet, like your validation_<tag>.csv:
#   one row per logger-hour, every main model's prediction as a column, plus
#   the inputs that produced it. Kept for every run in output/logs/.
validation <- preds %>% filter(model %in% c("M0", "M1", "M1d", "M2s", "M2r", "M3s", "M3r", "M4s", "M4r", "M4r_siteKd")) %>%
  select(model, sitekey, hoboID_env, datetime_pos, pred) %>%
  tidyr::pivot_wider(names_from = model, values_from = pred, names_prefix = "pred_") %>%
  left_join(eval_d %>% select(sitekey, hoboID_env, datetime_pos, light_mol, surface_par,
                              T_surf, water_depth_m, sun_altitude, kd_const, kd_sat,
                              kd_ref, shadow_refr, ambient),
            by = c("sitekey", "hoboID_env", "datetime_pos")) %>%
  left_join(preds %>% filter(model == "M4r") %>%
              select(sitekey, hoboID_env, datetime_pos, w_ambient_M4r = w_fold),
            by = c("sitekey", "hoboID_env", "datetime_pos")) %>%
  mutate(run_label = RUN_LABEL, run_id = RUN_ID,
         datetime_pos = format(datetime_pos, "%Y-%m-%d %H:%M:%S"))
val_file <- file.path(BASE_OUT_DIR, "logs", sprintf("validation_%s_%s.csv", RUN_ID, RUN_LABEL))
dir.create(dirname(val_file), showWarnings = FALSE, recursive = TRUE)
write.csv(validation, val_file, row.names = FALSE)
message("Validation sheet: ", val_file)

MAIN <- c("M0", "M1", "M1d", "M2s", "M2r", "M3s", "M3r", "M4s", "M4r")
STEPS <- list(c("M0", "M1"), c("M1", "M2s"), c("M1", "M2r"),
              c("M2s", "M3s"), c("M2r", "M3r"), c("M3s", "M4s"), c("M3r", "M4r"),
              c("M4r_siteKd", "M4r"),
              c("M1", "M1d"), c("M1d", "M3r"), c("M1d", "M4r"))   # [NEW] baseline steps
if (RUN_UNREFRACTED_DIAGNOSTIC) STEPS <- c(STEPS, list(c("M3r_unrefracted", "M3r")))

# ── 5. Metrics ───────────────────────────────────────────────────────────────
# [YOURS, EDITED] hourly and daytime-integrated metrics, per event and pooled
#   (your metrics_<RUN_TAG>.csv). [NEW] within-site metrics, events improved,
#   bootstrap CIs, ambient weight per fold
m_pooled <- preds %>% group_by(model) %>%
  group_modify(~ metrics(.x$light_mol, .x$pred)) %>% ungroup()
m_event  <- preds %>% group_by(model, sitekey) %>%
  group_modify(~ metrics(.x$light_mol, .x$pred)) %>% ungroup()
m_daytime <- preds %>% group_by(model) %>%
  group_modify(~ { di <- daytime_integrated(.x); metrics(di$obs, di$pred) }) %>% ungroup()
m_within <- preds %>% group_by(model) %>%
  group_modify(~ data.frame(within_r2 = within_r2(.x),
                            median_rank_rho = stats::median(within_rank(.x)$rho, na.rm = TRUE))) %>%
  ungroup()
w_folds <- preds %>% filter(!is.na(w_fold)) %>% distinct(model, sitekey, w_fold)

improved <- bind_rows(lapply(STEPS, function(s)
  data.frame(from = s[1], to = s[2], events_improved = events_improved(preds, s[1], s[2]),
             of = length(unique(preds$sitekey)))))

message("\nBootstrap (", BOOT_B, " replicates, loggers resampled within events)...")
boot <- bind_rows(lapply(intersect(MAIN, unique(preds$model)), function(m) {
  d <- preds[preds$model == m, ]
  fns <- list(rmse      = function(x) sqrt(mean((x$pred - x$light_mol)^2)),
              bias      = function(x) mean(x$pred - x$light_mol),
              within_r2 = within_r2)
  bind_rows(lapply(names(fns), function(nm)
    data.frame(model = m, metric = nm, t(boot_ci(d, fns[[nm]])))))
}))
# Note: the ambient weight is held at its cross-validated fold values inside the
# bootstrap (not re-fitted per replicate). Say so in the Methods.

write.csv(m_pooled,  out_path("metrics_pooled_hourly.csv"),  row.names = FALSE)
write.csv(m_event,   out_path("metrics_per_event.csv"),      row.names = FALSE)
write.csv(m_daytime, out_path("metrics_daytime_integrated.csv"), row.names = FALSE)
write.csv(m_within,  out_path("metrics_within_site.csv"),    row.names = FALSE)
write.csv(improved,  out_path("metrics_events_improved.csv"), row.names = FALSE)
write.csv(boot,      out_path("metrics_bootstrap_ci.csv"),   row.names = FALSE)
write.csv(w_folds,   out_path("ambient_weight_by_fold.csv"), row.names = FALSE)

message("\nPooled hourly metrics:")
print(as.data.frame(m_pooled), digits = 3, row.names = FALSE)
message("\nWithin-site structure:")
print(as.data.frame(m_within), digits = 3, row.names = FALSE)

# ── 6. How many reference loggers? ───────────────────────────────────────────
# [NEW]
message("\nReference-logger number experiment (", REF_N_DRAW, " draws per k)...")
exp_keys <- make_experiment_eval(eval_d, pool)
refn <- ref_number_experiment(cache, pool, exp_keys)
write.csv(refn, out_path("ref_number.csv"), row.names = FALSE)
print(refn %>% group_by(mode, k) %>%
        summarise(median_rmse = stats::median(rmse),
                  p05 = stats::quantile(rmse, 0.05), p95 = stats::quantile(rmse, 0.95),
                  .groups = "drop") %>% as.data.frame(), digits = 3)

pick <- function(tab, m, col) { v <- tab[[col]][tab$model == m]; if (length(v)) signif(v, 4) else NA }
headline <- c(M1d_rmse = pick(m_pooled, "M1d", "rmse"), M1d_r2 = pick(m_pooled, "M1d", "r2"),
              M4r_rmse = pick(m_pooled, "M4r", "rmse"), M4r_bias = pick(m_pooled, "M4r", "bias"),
              M4r_r2 = pick(m_pooled, "M4r", "r2"), M4r_within_r2 = pick(m_within, "M4r", "within_r2"),
              M3r_rmse = pick(m_pooled, "M3r", "rmse"), M3r_within_r2 = pick(m_within, "M3r", "within_r2"),
              M4r_dti_rmse = pick(m_daytime, "M4r", "rmse"),
              w_folds = paste(sort(unique(w_folds$w_fold[w_folds$model == "M4r"])), collapse = "/"))
write_manifest("03_models", extra = structure(c(
  sprintf("reference pool: %d loggers, %d used for Kd", nrow(pool), nrow(refs)),
  sprintf("evaluation rows: %d (frozen in eval_keys.csv)", N_EVAL),
  sprintf("ambient weight chosen per fold: %s",
          paste(sort(unique(w_folds$w_fold[w_folds$model == "M4r"])), collapse = ", ")),
  sprintf("run label: %s   ambient weight: %s", RUN_LABEL,
          if (is.finite(W_AMBIENT_FIXED)) paste("fixed", W_AMBIENT_FIXED) else "CV"),
  sprintf("validation sheet: %s", val_file)), headline = headline))
message("03 done.")
