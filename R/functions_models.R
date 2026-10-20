# ══════════════════════════════════════════════════════════════════════════════
#  functions_models.R — reference Kd, model ladder, cross-validation, metrics
#  Every model is arithmetic on the cached ray-tracing table, so no model
#  comparison ever needs new ray-tracing.
#
#  PROVENANCE TAGS (every function has one):
#    [YOURS]          your code from master_script.R, logic unchanged
#    [YOURS, EDITED]  your code, with the listed changes
#    [NEW]            written by Claude
#  Comments with no "[Claude]" marker are yours, copied from master_script.R.
#
#  [Claude] Most of this file is new. Removed from master_script.R:
#  w_ambient_combined() and its altitude/depth weighting, fit_alt_weights(),
#  fit_combined_weights(), fit_ambient_weights(), the per-site best-W block,
#  kd_implied, val_strat, find_reference_loggers() and the GT19 azimuth plot.
# ══════════════════════════════════════════════════════════════════════════════


# ── fresnel_T ─────────────────────────────────────────────────────────────────
# [NEW]
#   [Claude] Air-water transmittance for unpolarised light. Optional
#   (USE_FRESNEL in config.R); applied to predictions AND the reference-Kd
#   inversion together, otherwise Kd absorbs it.
fresnel_T <- function(alt_rad, n = N_WATER) {
  ti  <- pi / 2 - alt_rad
  out <- rep(1 - ((n - 1) / (n + 1))^2, length(ti))
  ok  <- ti > 1e-6
  tt  <- asin(sin(ti[ok]) / n)
  rs  <- (sin(ti[ok] - tt) / sin(ti[ok] + tt))^2
  rp  <- (tan(ti[ok] - tt) / tan(ti[ok] + tt))^2
  out[ok] <- 1 - (rs + rp) / 2
  out
}


# ── select_reference_pool ─────────────────────────────────────────────────────
# [NEW] — replaces your find_reference_loggers()
#   [Claude] Keeps your hard gates (min_hours, >= 2 midday hours) and your idea
#   of ranking by openness to the sun and sky. Drops the midday pred/meas
#   agreement term (it selected on the outcome) and replaces depth_sd < 3 m
#   (which every logger passed) with a minimum water depth. No measured light
#   is used, so selection cannot depend on how well the model fits.
select_reference_pool <- function(cache) {
  cache %>%
    mutate(hour = as.integer(format(datetime_pos, "%H"))) %>%
    group_by(sitekey, hoboID_env) %>%
    summarise(n_hours    = dplyr::n_distinct(datetime_pos),
              n_midday   = sum(hour %in% 10:14),
              shadow_min = min(shadow_refr, na.rm = TRUE),   # worst hour of the day
              ambient    = mean(ambient,    na.rm = TRUE),
              depth_min  = min(water_depth_m, na.rm = TRUE),
              .groups = "drop") %>%
    filter(n_hours   >= REF_MIN_HOURS,
           n_midday  >= REF_MIN_MIDDAY,
           ambient   >= REF_AMBIENT_MIN,
           shadow_min >= REF_SHADOW_MIN,
           is.finite(depth_min),                # [Claude] a logger with no depth can't give Kd
           depth_min >= REF_DEPTH_MIN_M) %>%
    group_by(sitekey) %>%
    arrange(desc(ambient), desc(shadow_min), .by_group = TRUE) %>%
    mutate(rank = dplyr::row_number()) %>%
    filter(rank <= REF_POOL_MAX) %>%
    ungroup()
}


# ── reference_kd ──────────────────────────────────────────────────────────────
# [YOURS, EDITED] — your derive_reference_kd()
#   - [Claude] works on the cached hourly table instead of re-averaging env_long,
#     so reference Kd and predictions use identical hourly values
#   - [Claude] hours with no valid value get the site's median REFERENCE Kd;
#     your run_pipeline() fell back to the satellite value for those hours
#   - [Claude] optional Fresnel factor (T_surf)
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
# if all 3 share the same bias.
#
# WHY MEDIAN NOT MEAN: median ignores one bad ref. Mean would drag the
# aggregate toward the bad ref.
reference_kd <- function(cache, ref_ids, quiet = FALSE) {
  r <- cache %>%
    semi_join(ref_ids, by = c("sitekey", "hoboID_env")) %>%
    mutate(transmit = light_mol / (surface_par * T_surf),
           kd_i = ifelse(is.finite(transmit) & transmit > 0 & transmit < 1 &
                           water_depth_m > 0,
                         -log(transmit) / water_depth_m, NA_real_))
  if (!quiet) {
    n_bad <- sum(is.na(r$kd_i))
    # transmit > 1 or wildly variable Kd means the ref logger isn't actually
    # open / has a unit mismatch with surface_par / has a wrong depth
    if (n_bad > 0) message("  [ref-Kd] ", n_bad, " of ", nrow(r),
                           " reference logger-hours invalid (transmittance outside 0-1)")
  }
  hourly <- r %>% group_by(sitekey, datetime_pos) %>%
    summarise(kd_ref_hour = if (all(is.na(kd_i))) NA_real_ else stats::median(kd_i, na.rm = TRUE),
              n_ref = sum(is.finite(kd_i)),   # how many refs actually contributed each hour
              .groups = "drop")
  site_med <- hourly %>% group_by(sitekey) %>%
    summarise(kd_ref_site = stats::median(kd_ref_hour, na.rm = TRUE), .groups = "drop")
  cache %>%
    distinct(sitekey, datetime_pos) %>%
    left_join(hourly,   by = c("sitekey", "datetime_pos")) %>%
    left_join(site_med, by = "sitekey") %>%
    mutate(kd_filled = !is.finite(kd_ref_hour),
           kd_ref    = ifelse(kd_filled, kd_ref_site, kd_ref_hour))
}


# ── spectral_atten ────────────────────────────────────────────────────────────
# [YOURS, EDITED] — your calc_water_attenuation_spectral()
#   - [Claude] vectorised over sites; for the supplementary model only
#
# Depth-dependent spectral path. Instead of one KdPAR, attenuate each band at its
# own Kd and sum the surviving PAR fraction at depth; depth-dependence emerges
# automatically (Canto method). f_lambda = surface photon weight per band (∝ λ,
# same assumption as kdpar_photon; Kirk 2011). Drops NA bands and renormalises so
# a site missing one band (e.g. CB18, 9 bands) still computes.
spectral_atten <- function(depth, sitekey, kd_bands) {
  out <- rep(NA_real_, length(depth))
  for (k in unique(sitekey)) {
    b <- kd_bands[kd_bands$sitekey == k & is.finite(kd_bands$kd), ]
    if (nrow(b) < 2) next   # need at least 2 bands for a meaningful sum
    f <- b$wave_nm / sum(b$wave_nm)
    i <- which(sitekey == k)
    out[i] <- vapply(depth[i], function(z)
      if (!is.finite(z) || z <= 0) NA_real_ else sum(f * exp(-b$kd * z)), numeric(1))
  }
  out
}


# ── predict_par ───────────────────────────────────────────────────────────────
# [NEW] — but the equation is yours, from process_timestep():
#      pred_par = surface_par × (w_direct·shadow + w_ambient·ambient) × water_atten
#   [Claude] One function for every rung of the model ladder, with one static
#   weight. kd_col = NULL means no attenuation (M0).
#
# WHY w_direct + w_ambient = 1:
# These are the PARTITION of incoming surface PAR into its direct-beam and
# diffuse-sky components — physically exhaustive (no third component), so the
# two fractions sum to 1 by definition. This keeps eff_light in [0,1] and gives
# it a clean meaning: "fraction of surface PAR surviving geometry (beam-block +
# sky-occlusion) before depth attenuation."
predict_par <- function(d, kd_col = NULL,
                        shading = c("none", "direct", "direct_air", "both"),
                        w = 0, shadow_col = "shadow_refr", ambient_col = "ambient") {
  shading <- match.arg(shading)
  eff <- switch(shading,
                none       = 1,
                direct     = d[[shadow_col]],
                direct_air = d$shadow_air,
                both       = (1 - w) * d[[shadow_col]] + w * d[[ambient_col]])
  # Beer-Lambert law: exp(−Kd × depth)
  att <- if (is.null(kd_col)) 1 else exp(-d[[kd_col]] * d$water_depth_m)
  d$surface_par * d$T_surf * eff * att
}


# ── loeo_w ────────────────────────────────────────────────────────────────────
# [NEW] — the grid is the same as your per-site best-W block (0 to 0.95 by 0.05;
#   here 0 to 1), but the weight is chosen on 8 events and scored on the 9th
#   [Claude] Returns out-of-fold predictions (`pred`) and each fold's weight.
loeo_w <- function(d, kd_col, shadow_col = "shadow_refr", ambient_col = "ambient",
                   grid = W_GRID) {
  bind_rows(lapply(unique(d$sitekey), function(k) {
    train <- d[d$sitekey != k, ]
    test  <- d[d$sitekey == k, ]
    rmse  <- vapply(grid, function(w)
      sqrt(mean((predict_par(train, kd_col, "both", w, shadow_col, ambient_col) -
                   train$light_mol)^2)), numeric(1))
    w_best <- grid[which.min(rmse)]
    test$pred   <- predict_par(test, kd_col, "both", w_best, shadow_col, ambient_col)
    test$w_fold <- w_best
    test
  }))
}


# ── metrics ───────────────────────────────────────────────────────────────────
# [YOURS] — your metrics helper
#   ([Claude] r2 written as cor()^2, which equals your summary(lm)$r.squared)
metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  n <- length(obs)
  if (n < 3) return(data.frame(n = n, rmse = NA, mae = NA, bias = NA,
                               r2 = NA, r2_1to1 = NA, slope = NA, intercept = NA))
  res <- pred - obs
  fit <- stats::lm(pred ~ obs)
  data.frame(n = n,
             rmse = sqrt(mean(res^2)),
             mae  = mean(abs(res)),
             bias = mean(res),                                   # +ve = over-predicts
             r2   = cor(obs, pred)^2,                            # Pearson R² (best-fit line)
             r2_1to1 = 1 - sum(res^2) / sum((obs - mean(obs))^2), # Nash-Sutcliffe vs 1:1
             slope = unname(stats::coef(fit)[2]),
             intercept = unname(stats::coef(fit)[1]))
}


# ── daytime_integrated ────────────────────────────────────────────────────────
# [YOURS, EDITED] — your DLI block
#   - [Claude] renamed: only hours with sun altitude above SUN_ALT_MIN_DEG are
#     included, so it is daytime-integrated PAR, not DLI
#
#  integrate BOTH from the same table, which holds ONE hourly-mean
#  row per logger/hour for both light_mol and pred_par. Both are mean
#  PAR over the hour in umol m^-2 s^-1, so BOTH multiply by 3600 s and divide by
#  1e6 -> mol m^-2. Identical rows in, identical window out: directly comparable.
daytime_integrated <- function(d, pred_col = "pred") {
  d %>% group_by(sitekey, hoboID_env) %>%
    summarise(obs   = sum(light_mol)        * 3600 / 1e6,
              pred  = sum(.data[[pred_col]]) * 3600 / 1e6,
              depth = mean(water_depth_m),
              .groups = "drop")
}


# ── within_r2 / within_anomalies / within_rank ────────────────────────────────
# [NEW]
#   [Claude] Remove the site-hour mean from measured and predicted PAR, then
#   score what is left. Only logger-to-logger differences count, which is what
#   reef structure is supposed to explain.
within_r2 <- function(d, pred_col = "pred", min_loggers = 3) {
  a <- within_anomalies(d, pred_col, min_loggers)
  1 - sum((a$d_obs - a$d_pred)^2) / sum(a$d_obs^2)
}

within_anomalies <- function(d, pred_col = "pred", min_loggers = 3) {
  d %>% group_by(sitekey, datetime_pos) %>%
    filter(dplyr::n() >= min_loggers) %>%
    mutate(d_obs  = light_mol - mean(light_mol),
           d_pred = .data[[pred_col]] - mean(.data[[pred_col]])) %>%
    ungroup()
}

within_rank <- function(d, pred_col = "pred") {
  daytime_integrated(d, pred_col) %>%
    group_by(sitekey) %>%
    summarise(rho = suppressWarnings(stats::cor(obs, pred, method = "spearman")),
              .groups = "drop")
}


# ── boot_ci ───────────────────────────────────────────────────────────────────
# [NEW]
#   [Claude] Cluster bootstrap: resample loggers within each sampling event.
#   `fun` takes a data frame and returns ONE number. Duplicated loggers get a
#   unique id so within-site statistics treat them as separate loggers.
boot_ci <- function(d, fun, B = BOOT_B, seed = SEED) {
  set.seed(seed)
  loggers <- distinct(d, sitekey, hoboID_env)
  vals <- replicate(B, {
    s <- loggers %>% group_by(sitekey) %>%
      slice_sample(prop = 1, replace = TRUE) %>% ungroup() %>%
      mutate(boot_id = paste0(hoboID_env, "_", dplyr::row_number()))
    bd <- inner_join(s, d, by = c("sitekey", "hoboID_env"),
                     relationship = "many-to-many") %>%
      mutate(hoboID_env = boot_id)
    fun(bd)
  })
  stats::setNames(stats::quantile(vals, c(0.025, 0.5, 0.975), na.rm = TRUE),
                  c("lo", "mid", "hi"))
}


# ── events_improved ───────────────────────────────────────────────────────────
# [NEW]
#   [Claude] In how many sampling events does model b have lower RMSE than a.
#   Replaces the paired t-tests / Wilcoxon on pooled hours.
events_improved <- function(preds, a, b) {
  per <- preds %>% filter(model %in% c(a, b)) %>%
    group_by(model, sitekey) %>%
    summarise(rmse = sqrt(mean((pred - light_mol)^2)), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = model, values_from = rmse)
  sum(per[[b]] < per[[a]], na.rm = TRUE)   # [Claude] events scored by both models
}


# ── make_experiment_eval / ref_number_experiment ──────────────────────────────
# [NEW]
#   [Claude] How many reference loggers are needed, and does choosing open
#   positions matter? Holds back a seeded random fraction of non-pool loggers
#   for evaluation, so "random" draws have loggers to choose from that are not
#   being scored.
make_experiment_eval <- function(cache, pool, frac = 0.5, seed = SEED) {
  set.seed(seed)
  keep <- cache %>% distinct(sitekey, hoboID_env) %>%
    anti_join(pool, by = c("sitekey", "hoboID_env")) %>%
    group_by(sitekey) %>% slice_sample(prop = frac) %>% ungroup()
  cache %>% semi_join(keep, by = c("sitekey", "hoboID_env")) %>%
    distinct(sitekey, hoboID_env, datetime_pos)
}

ref_number_experiment <- function(cache, pool, eval_keys,
                                  ks = REF_KS, n_draw = REF_N_DRAW, seed = SEED) {
  set.seed(seed)
  eval_d <- semi_join(cache, eval_keys, by = c("sitekey", "hoboID_env", "datetime_pos")) %>%
    select(-any_of(c("kd_ref", "kd_ref_site", "kd_filled")))   # replaced per draw below
  non_eval <- cache %>% distinct(sitekey, hoboID_env) %>%
    anti_join(distinct(eval_keys, sitekey, hoboID_env), by = c("sitekey", "hoboID_env"))
  out <- list()
  for (mode in c("open", "random")) for (k in ks) for (i in seq_len(n_draw)) {
    src <- if (mode == "open") pool else non_eval
    ids <- src %>% group_by(sitekey) %>%
      slice_sample(n = k) %>% ungroup() %>% select(sitekey, hoboID_env)
    kd  <- reference_kd(cache, ids, quiet = TRUE)
    dd  <- left_join(eval_d, select(kd, sitekey, datetime_pos, kd_ref),
                     by = c("sitekey", "datetime_pos"))
    p   <- predict_par(dd, "kd_ref", "direct")
    ok  <- is.finite(p)
    out[[length(out) + 1]] <- data.frame(
      mode = mode, k = k, draw = i,
      rmse = sqrt(mean((p[ok] - dd$light_mol[ok])^2)),
      n_eval = sum(ok))
  }
  bind_rows(out)
}
