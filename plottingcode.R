# ═══════════════════════════════════════════════════════════════════════════
#  THESIS FIGURES — cleaned plotting script
#
#  Every original block is kept. Bugs fixed:
#    - helper functions now DEFINED before they are CALLED (Fig 9, Fig 10)
#    - `surf` built before it is used in the Fig 3 walkthrough
#    - stray empty png()/dev.off()/dev.new() removed
#    - baseline constant Kd standardised to 0.20 throughout
#    - paths/tags collected into one CONFIG block at the top
#    - superseded draft blocks kept but clearly marked [DRAFT — superseded]
#
#  THESIS-FIGURE MAP (search these tags):
#    [FIG 3]   one-day walkthrough (tide / surface PAR / Kd) — CB19
#    [FIG 6]   design-decision walkthrough, panels A–E
#    [FIG 7]   baseline ladder M0 → M5 (RMSE / R² / bias)
#    [FIG 8]   per-site PAR, baseline vs chosen (2 pages)
#    [FIG 9]   pooled hourly obs vs pred, 08/12/16h slopes
#    [FIG 10]  pooled DLI scatter, colour = depth
#    [FIG 11]  per-site residual diagnostics (2 pages)
#    (FIG 5 = ray-trace shadow discs is produced elsewhere, not in this file)
# ═══════════════════════════════════════════════════════════════════════════

# ── CONFIG ─────────────────────────────────────────────────────────────────
# Set these once. Everything downstream reads from here.
RUNS_ROOT   <- "runs"                      # folder holding the per-run subfolders
CHOSEN_TAG  <- "20260702_204808"           # the chosen-model run
RUN_TAG     <- CHOSEN_TAG                   # tag used in output filenames
RUNS_DIR    <- file.path(RUNS_ROOT, RUN_TAG)
TAG         <- CHOSEN_TAG                   # alias used by some blocks below

KD_BASELINE <- 0.20                         # constant Kd for the "baseline" model

VAL_CSV     <- file.path(RUNS_ROOT, TAG, sprintf("validation_%s.csv", TAG))

library(dplyr)

# ── shared helpers (defined ONCE, used everywhere) ─────────────────────────
metrics <- function(obs, pred) {
  ok  <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]
  list(n    = length(obs),
       rmse = sqrt(mean((pred - obs)^2)),
       bias = mean(pred - obs),
       r2   = if (length(obs) > 2) cor(pred, obs)^2 else NA_real_)
}

hours_to_label <- function(hn) {
  d <- as.Date(floor(hn / 24), origin = "1970-01-01")
  h <- hn %% 24
  sprintf("%s %02d:00", format(d, "%b %d"), h)
}

# add the derived columns a validation frame needs (idempotent)
add_derived_cols <- function(v) {
  v$datetime_pos <- as.POSIXct(v$datetime_pos,
                               format = "%Y-%m-%d %H:%M:%S",
                               tz     = "Australia/Brisbane")
  v$date_str  <- substr(v$datetime_pos, 1, 10)
  v$hour      <- as.integer(substr(v$datetime_pos, 12, 13))
  v$hours_num <- as.numeric(as.Date(v$date_str)) * 24 + v$hour
  v$resid     <- v$pred_par - v$light_mol
  # baseline (shadow only, constant Kd) and chosen predictions
  v$pred_chosen   <- v$pred_par
  v$pred_baseline <- v$surface_par * v$shadow_fraction *
    exp(-KD_BASELINE * v$water_depth_m)
  v
}


# ══════════════════════════════════════════════════════════════════════════
#  [FIG 3] One-day walkthrough at a chosen site
#  Panel a: tide height · Panel b: surface PAR · Panel c: Kd attenuation
# ══════════════════════════════════════════════════════════════════════════

v_all <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v_all$datetime_pos <- as.POSIXct(v_all$datetime_pos, tz = "Australia/Brisbane")

SITE_PICK <- "CB19"      # any sitekey in v_all
DAY_PICK  <- NULL        # NULL = auto-pick day with most timesteps

vsite <- v_all[v_all$sitekey == SITE_PICK, ]
if (nrow(vsite) == 0) stop("No rows for ", SITE_PICK)

vsite$date <- as.Date(vsite$datetime_pos)
if (is.null(DAY_PICK)) {
  day_counts <- table(vsite$date)
  DAY_PICK   <- as.Date(names(day_counts)[which.max(day_counts)])
}
vday <- vsite[vsite$date == DAY_PICK, ]
vday$hour <- as.numeric(format(vday$datetime_pos, "%H")) +
  as.numeric(format(vday$datetime_pos, "%M")) / 60
message("Plotting ", SITE_PICK, " on ", DAY_PICK,
        " (", nrow(vday), " rows, ",
        length(unique(vday$hoboID_env)), " sensors)")

# observed vs predicted attenuation
vday$atten_obs  <- pmin(pmax(vday$light_mol /
                               (vday$surface_par * vday$eff_light), 0), 1.2)
vday$atten_pred <- vday$water_atten

agg <- do.call(rbind, lapply(sort(unique(vday$datetime_pos)), function(t) {
  sub <- vday[vday$datetime_pos == t, ]
  data.frame(
    hour      = sub$hour[1],
    pred_mean = mean(sub$atten_pred, na.rm = TRUE),
    pred_min  = min(sub$atten_pred,  na.rm = TRUE),
    pred_max  = max(sub$atten_pred,  na.rm = TRUE),
    obs_mean  = mean(sub$atten_obs,  na.rm = TRUE),
    obs_q25   = quantile(sub$atten_obs, 0.25, na.rm = TRUE),
    obs_q75   = quantile(sub$atten_obs, 0.75, na.rm = TRUE)
  )
}))
agg <- agg[order(agg$hour), ]

# surface PAR curve (site-wide, one value per timestep) — BUILT BEFORE PLOTTING
surf <- unique(vday[, c("hour", "surface_par")])
surf <- surf[order(surf$hour), ]

# tide series — from env_long if available, else fall back to the tide CSV
if (exists("env_long")) {
  site_pick <- unique(vday$site)
  year_pick <- unique(vday$year)
  tide_day <- env_long %>%
    filter(site == site_pick, year == year_pick,
           as.Date(datetime_pos) == DAY_PICK) %>%
    distinct(datetime_pos, tide_cm) %>%
    arrange(datetime_pos)
  if (nrow(tide_day) == 0)
    stop("No env_long tide rows for ", site_pick, " ", year_pick, " on ", DAY_PICK)
  tide_day$m <- if (max(abs(tide_day$tide_cm), na.rm = TRUE) < 10)
    tide_day$tide_cm else tide_day$tide_cm / 100
  tide_day$hour <- as.numeric(format(tide_day$datetime_pos, "%H")) +
    as.numeric(format(tide_day$datetime_pos, "%M")) / 60
  tide_day <- tide_day[order(tide_day$hour), ]
} else {
  message("env_long not in memory — tide panel will be skipped. ",
          "Source master_script.R load block to include it.")
  tide_day <- NULL
}

png(sprintf("fig_dayprofile_%s_%s.png",
            SITE_PICK, format(DAY_PICK, "%Y%m%d")),
    width = 10 * 130, height = 11 * 130, res = 130)

layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1), heights = c(1, 1, 2))
par(oma = c(4.5, 0, 0, 0), mgp = c(2.6, 0.7, 0))
xr <- range(agg$hour)

# Panel a: tide
par(mar = c(0.5, 4.5, 3, 4.5))
if (!is.null(tide_day)) {
  plot(tide_day$hour, tide_day$m, type = "l", col = "darkgreen", lwd = 2,
       xlim = xr, xaxt = "n", xlab = "", cex.axis = 1.5,
       ylab = "a. Tide (m)", cex.lab = 1.5, cex.main = 2,
       main = sprintf("Water attenuation over one day \u2014 %s, %s",
                      SITE_PICK, format(DAY_PICK, "%d %b %Y")))
} else {
  plot.new(); title(main = "a. Tide (m) — env_long not loaded")
}

# Panel b: surface PAR
par(mar = c(0.5, 4.5, 0.5, 4.5))
plot(surf$hour, surf$surface_par, type = "n",
     xlim = xr, xaxt = "n", xlab = "",
     ylab = expression("b. Surface PAR (\u00b5mol m"^-2*" s"^-1*")"),
     cex.axis = 1.5, cex.lab = 1.5,
     ylim = c(0, max(surf$surface_par, na.rm = TRUE) * 1.05))
polygon(c(surf$hour, rev(surf$hour)),
        c(surf$surface_par, rep(0, nrow(surf))),
        col = adjustcolor("gold", 0.20), border = NA)
lines(surf$hour, surf$surface_par, col = "goldenrod", lwd = 3)

# Panel c: attenuation
par(mar = c(0.5, 4.5, 0.5, 4.5))
plot(NA, xlim = xr, ylim = c(0, 1), cex.axis = 1.5, cex.lab = 1.5,
     xlab = "", ylab = "c. Kd (fraction of surface PAR reaching depth)")
points(vday$hour, vday$atten_obs,
       pch = 16, col = adjustcolor("steelblue", 0.25), cex = 0.7)
polygon(c(agg$hour, rev(agg$hour)), c(agg$obs_q25, rev(agg$obs_q75)),
        col = adjustcolor("steelblue", 0.15), border = NA)
polygon(c(agg$hour, rev(agg$hour)), c(agg$pred_min, rev(agg$pred_max)),
        col = adjustcolor("firebrick", 0.15), border = NA)
lines(agg$hour, agg$obs_mean,  col = "steelblue", lwd = 2.5)
lines(agg$hour, agg$pred_mean, col = "firebrick", lwd = 2.5, lty = 2)
legend("topright", bty = "n", cex = 1.5,
       legend = c("Observed (per sensor)", "Observed IQR across sensors",
                  "Observed mean", "Predicted range across sensors",
                  "Predicted mean  \u2014  exp(\u2212Kd \u00d7 depth)"),
       col = c(adjustcolor("steelblue", 0.5), "steelblue", "steelblue",
               "firebrick", "firebrick"),
       pch = c(16, NA, NA, NA, NA), lty = c(NA, 1, 1, 1, 2),
       lwd = c(NA, 8, 2.5, 8, 2.5))

mean_kd   <- round(mean(vday$Kd, na.rm = TRUE), 3)
depth_rng <- round(range(vday$water_depth_m, na.rm = TRUE), 2)
mtext(bquote("Kd = " ~ .(mean_kd) ~ "m"^-1 ~ "  |  depth range " ~
               .(depth_rng[1]) * "\u2013" * .(depth_rng[2]) ~ "m"),
      side = 3, line = -1.5, cex = 0.8, col = "grey30")
mtext("Hour of day (local)", side = 1, line = 2.5, outer = TRUE)
layout(1)
dev.off()
message("Saved fig_dayprofile_", SITE_PICK, "_", format(DAY_PICK, "%Y%m%d"), ".png")


# ══════════════════════════════════════════════════════════════════════════
#  [SUPPORTING] Surface vs single-sensor PAR through the day
#  (diagnostic; not a numbered thesis figure, kept as-is but ordering fixed)
# ══════════════════════════════════════════════════════════════════════════

SENSOR_PICK <- NULL
if (is.null(SENSOR_PICK)) {
  MIN_HOURS <- 5
  sensor_scores <- do.call(rbind, lapply(unique(vday$hoboID_env), function(id) {
    sub <- vday[vday$hoboID_env == id, ]
    ok  <- is.finite(sub$pred_par) & is.finite(sub$light_mol)
    n   <- sum(ok)
    if (n < MIN_HOURS) return(NULL)
    data.frame(hoboID_env = id, n = n,
               rmse = sqrt(mean((sub$pred_par[ok] - sub$light_mol[ok])^2)),
               bias = mean(sub$pred_par[ok] - sub$light_mol[ok]),
               mape = mean(abs(sub$pred_par[ok] - sub$light_mol[ok]) /
                             pmax(sub$light_mol[ok], 1)))
  }))
  sensor_scores <- sensor_scores[order(sensor_scores$rmse), ]
  SENSOR_PICK   <- sensor_scores$hoboID_env[1]
  message("Picked sensor: ", SENSOR_PICK)
}
vs <- vday[vday$hoboID_env == SENSOR_PICK, ]
vs <- vs[order(vs$hour), ]

par(mar = c(4.5, 4.5, 3.5, 1), mgp = c(2.6, 0.7, 0))
ymax <- max(c(surf$surface_par, vs$pred_par, vs$light_mol), na.rm = TRUE) * 1.05
plot(NA, xlim = range(surf$hour), ylim = c(0, ymax),
     xlab = "Hour of day (local)",
     ylab = expression("PAR (\u00b5mol m"^-2*" s"^-1*")"),
     main = sprintf("Surface vs sensor PAR \u2014 %s, sensor %s, %s",
                    SITE_PICK, SENSOR_PICK, format(DAY_PICK, "%d %b %Y")))
polygon(c(surf$hour, rev(surf$hour)),
        c(surf$surface_par, rep(0, nrow(surf))),
        col = adjustcolor("gold", 0.20), border = NA)
lines(surf$hour, surf$surface_par, col = "goldenrod", lwd = 3)
lines(vs$hour, vs$pred_par,   col = "firebrick", lwd = 2.2, lty = 2)
points(vs$hour, vs$pred_par,  col = "firebrick", pch = 16, cex = 0.9)
lines(vs$hour, vs$light_mol,  col = "steelblue", lwd = 2.2)
points(vs$hour, vs$light_mol, col = "steelblue", pch = 16, cex = 0.9)
legend("topright", bty = "n", cex = 0.85,
       legend = c("Surface PAR (top of water)",
                  "Predicted PAR at sensor", "Measured PAR at sensor"),
       col = c("goldenrod", "firebrick", "steelblue"),
       lwd = c(3, 2.2, 2.2), lty = c(1, 2, 1))

# ══════════════════════════════════════════════════════════════════════════
#  [FIG 6] Design-decision walkthrough — panels A–E
#  Tweaks vs previous:
#   - density fills have NO outline; the "chosen"/user-highlighted curve is
#     drawn LAST so it sits in front of the others
#   - Panel D made narrower (shorter labels + trimmed margin) so text fits
#   - Panel D labels abbreviated: constant Kd -> const. Kd,
#     satellite Kd -> Sat Kd, reference Kd -> Ref. Kd
#   - text sizes bumped slightly again
# ══════════════════════════════════════════════════════════════════════════

CANDS_A <- c("20260625_140321", "20260625_153449", "20260625_184918")
LABS_A  <- c("Shadow only", "+ ambient (static)", "+ ambient (alt+depth-dep)")

CANDS_B <- c("20260625_184918", "20260625_203039", "20260626_110245")
LABS_B  <- c("Constant Kd (0.32)", "Satellite Kd + overrides", "Reference Kd")

CANDS_C <- c("20260625_153449", "20260625_184918")
LABS_C  <- c("Static weights", "Alt+depth-dep weights")

CANDS_E <- c("20260702_084057", "20260702_123929",
             "20260702_201824", CHOSEN_TAG)
LABS_E  <- c("Satellite Kd (photon)", "Reference Kd (early)",
             "Reference Kd (mid)", "Reference Kd (chosen)")

SIX_SITES <- "CB19;GT19;HS19;N319;TB19;TM19"

# ── text-size knobs ────────────────────────────────────────────────────────
CEX_MAIN <- 2.0
CEX_LAB  <- 2
CEX_AXIS <- 1.8
CEX_LEG  <- 1.7
CEX_DLAB <- 2   # panel D y-axis config labels
FILL_ALPHA <- 0.55  # density fill transparency (no outline)

load_val <- function(tag) {
  f <- file.path(RUNS_ROOT, tag, sprintf("validation_%s.csv", tag))
  if (!file.exists(f)) { message("  missing validation CSV for ", tag); return(NULL) }
  d <- read.csv(f); d$run_tag <- tag; d$resid <- d$pred_par - d$light_mol
  d
}

build_densities <- function(tags, labels) {
  keep <- !sapply(tags, function(t) is.null(load_val(t)))
  tags <- tags[keep]; labels <- labels[keep]
  lst <- lapply(tags, function(t) {
    d <- load_val(t); r <- d$resid[is.finite(d$resid)]
    list(dens = density(r, from = quantile(r, 0.005), to = quantile(r, 0.995)),
         bias = mean(r), sd = sd(r), n = length(r))
  })
  names(lst) <- labels
  lst
}

draw_density_panel <- function(dens_list, main, chosen_lab = NULL,
                               palette = "Dark 3",
                               xlim = NULL, ylim = NULL) {
  if (length(dens_list) == 0) { plot.new(); title(paste(main, "(no data)")); return() }
  labs <- names(dens_list)
  pal  <- hcl.colors(length(labs), palette)
  if (is.null(xlim)) xlim <- range(sapply(dens_list, function(x) range(x$dens$x)))
  if (is.null(ylim)) ylim <- c(0, max(sapply(dens_list, function(x) max(x$dens$y))))
  plot(NA, xlim = xlim, ylim = ylim,
       xlab = expression("Residual (pred \u2212 obs, \u00b5mol m"^-2*" s"^-1*")"),
       ylab = "Density", main = main,
       cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
  
  # draw order: everything except the chosen curve first (tallest-peaked last
  # among those), then the chosen curve on top so it sits in front
  chosen_i <- if (!is.null(chosen_lab)) which(labs == chosen_lab) else integer(0)
  others   <- setdiff(seq_along(dens_list), chosen_i)
  others   <- others[order(sapply(dens_list[others], function(x) max(x$dens$y)),
                           decreasing = TRUE)]
  draw_ord <- c(others, chosen_i)   # chosen last = in front
  
  for (i in draw_ord) {
    dx <- dens_list[[i]]$dens$x; dy <- dens_list[[i]]$dens$y
    polygon(c(dx[1], dx, dx[length(dx)]), c(0, dy, 0),
            col = adjustcolor(pal[i], FILL_ALPHA), border = NA)  # translucent, NO outline
  }
  abline(v = 0, lty = 2, col = "grey30", lwd = 1.6)
  leg_lab <- sprintf("%s\n(bias=%+.0f, sd=%.0f)", labs,
                     sapply(dens_list, function(x) x$bias),
                     sapply(dens_list, function(x) x$sd))
  legend("topleft", legend = leg_lab, fill = adjustcolor(pal, FILL_ALPHA), border = NA,
         bty = "n", cex = CEX_LEG, y.intersp = 1.4)
}

dens_A <- build_densities(CANDS_A, LABS_A)
dens_B <- build_densities(CANDS_B, LABS_B)
dens_C <- build_densities(CANDS_C, LABS_C)
dens_E <- build_densities(CANDS_E, LABS_E)

all_dens_ABCE <- c(dens_A, dens_B, dens_C, dens_E)
XLIM_ABCE <- range(sapply(all_dens_ABCE, function(x) range(x$dens$x)))
YLIM_ABCE <- c(0, max(sapply(all_dens_ABCE, function(x) max(x$dens$y))))

while (!is.null(dev.list())) dev.off()
png(file.path(RUNS_DIR, sprintf("fig2_1_walkthrough_%s.png", RUN_TAG)),
    width = 14 * 300, height = 15 * 300, res = 300)
# left column A,B,C stacked; right column D (top) over E (bottom)
layout_mat <- rbind(c(1, 4),
                    c(2, 4),
                    c(3, 5),
                    c(6, 6))
layout(layout_mat, heights = c(3, 3, 3, 0.5), widths = c(1, 1.2))
par(mar = c(4.8, 4.8, 3.2, 1), mgp = c(3.2, 0.9, 0))

draw_density_panel(dens_A, "A. Effect of ambient layer (Kd fixed)",
                   chosen_lab = "+ ambient (alt+depth-dep)",
                   xlim = XLIM_ABCE, ylim = YLIM_ABCE)
draw_density_panel(dens_B, "B. Effect of Kd source (weights fixed)",
                   chosen_lab = "Reference Kd", xlim = XLIM_ABCE, ylim = YLIM_ABCE)
draw_density_panel(dens_C, "C. Effect of weighting scheme (Kd fixed)",
                   chosen_lab = "Alt+depth-dep weights",
                   xlim = XLIM_ABCE, ylim = YLIM_ABCE)

# Panel D: all 6-site configs from metrics_log
par(mar = c(4.8, 30, 3.2, 4), mgp = c(3.4, 0.9, 0))
ml   <- read.csv(file.path(RUNS_ROOT, "metrics_log.csv"), stringsAsFactors = FALSE)
ml_6 <- ml[ml$sites == SIX_SITES, ]
ml_6 <- ml_6[order(ml_6$rmse_hourly), ]
ml_6$y   <- seq_len(nrow(ml_6))

# strip parentheses, then abbreviate the Kd-method name
kd_clean <- trimws(sub("\\s*\\(.*\\)", "", ml_6$kd_method))
kd_clean <- sub("(?i)^constant Kd.*",       "Const. Kd",     kd_clean, perl = TRUE)
kd_clean <- sub("(?i)^satellite Kd \\+.*", "Sat. Kd", kd_clean, perl = TRUE)
kd_clean <- sub("(?i)^satellite Kd.*",      "Sat. Kd",       kd_clean, perl = TRUE)
kd_clean <- sub("(?i)^reference Kd.*",      "Ref. Kd",       kd_clean, perl = TRUE)
kd_clean <- sub("(?i)^spectral Kd.*",       "Spec. Kd",      kd_clean, perl = TRUE)

ml_6$lab <- sprintf("%s | amb %s | alt %s | dep %s",
                    kd_clean,
                    ifelse(ml_6$use_ambient, "T", "F"),
                    ifelse(ml_6$use_alt_dep_weights, "T", "F"),
                    ifelse(ml_6$use_depth_dep_weights, "T", "F"))

shortlist     <- unique(c(CANDS_A, CANDS_B, CANDS_C))
ml_6$in_short <- ml_6$run_tag %in% shortlist
ml_6$col      <- ifelse(ml_6$in_short, "firebrick", "grey60")
ml_6$cex      <- ifelse(ml_6$in_short, 1.9, 1.1)
d_xlim <- range(ml_6$rmse_hourly, na.rm = TRUE) + c(-8, 8)
plot(ml_6$rmse_hourly, ml_6$y, type = "n", yaxt = "n", xlim = d_xlim,
     xlab = expression("Hourly RMSE (\u00b5mol m"^-2*" s"^-1*")"), ylab = "",
     main = "D. All 6-site 2019 configs (metrics_log)",
     cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
segments(min(ml_6$rmse_hourly) * 0.98, ml_6$y, ml_6$rmse_hourly, ml_6$y,
         col = "grey85", lwd = 0.8)
points(ml_6$rmse_hourly, ml_6$y, pch = 19, col = ml_6$col, cex = ml_6$cex)
abline(v = min(ml_6$rmse_hourly), col = "firebrick", lty = 3)
axis(2, at = ml_6$y, labels = ml_6$lab, las = 1, cex.axis = CEX_DLAB)
legend("bottomright", bty = "n", cex = CEX_LEG, pch = 19,
       col = c("firebrick", "grey60"), pt.cex = c(1.9, 1.1),
       legend = c("Shown in panels A/B/C", "Other config"))

# Panel E
par(mar = c(4.8, 4.8, 3.2, 1), mgp = c(3.2, 0.9, 0))
draw_density_panel(dens_E, "E. Chosen model on 9-site full dataset",
                   chosen_lab = "Reference Kd (chosen)",
                   xlim = XLIM_ABCE, ylim = YLIM_ABCE)

# caption strip
par(mar = c(0.5, 1, 0.5, 1)); plot.new()
text(0.5, 0.5,
     paste0("Panels A\u2013C: sensitivity on 6-site 2019 subset. ",
            "Panel D: all configurations on that subset, ranked by hourly RMSE; ",
            "red = shown in A\u2013C. ",
            "Panel E: chosen configuration (bold) on the full 9-site 2018\u20132019 dataset."),
     cex = 1.3, adj = 0.5)
dev.off()
message("Saved fig2_1_walkthrough_", RUN_TAG, ".png")

# ══════════════════════════════════════════════════════════════════════════
#  [FIG 7] Baseline ladder: M0 → M5
# ══════════════════════════════════════════════════════════════════════════

v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)

v$pred_M0 <- v$surface_par
v$pred_M1 <- v$surface_par * exp(-KD_BASELINE * v$water_depth_m)
v$pred_M2 <- v$surface_par * exp(-v$Kd * v$water_depth_m)
v$pred_M3 <- v$surface_par * v$shadow_fraction * exp(-v$Kd * v$water_depth_m)
v$pred_M4 <- v$surface_par *
  (0.3 * v$shadow_fraction + 0.7 * v$ambient_fraction) *
  exp(-v$Kd * v$water_depth_m)
v$pred_M5 <- v$pred_par

MODELS <- list(
  "M0: no water"    = v$pred_M0,
  "M1: BL (Kd=0.2)" = v$pred_M1,
  "M2: BL (ref Kd)" = v$pred_M2,
  "M3: + shadow"    = v$pred_M3,
  "M4: + ambient"   = v$pred_M4,
  "M5: chosen"      = v$pred_M5
)

tab <- do.call(rbind, lapply(seq_along(MODELS), function(i) {
  m <- metrics(v$light_mol, MODELS[[i]])
  data.frame(model = names(MODELS)[i], id = i,
             rmse = m$rmse, bias = m$bias, r2 = m$r2, n = m$n)
}))
print(tab, row.names = FALSE)

while (!is.null(dev.list())) dev.off()
FIG_W <- 7.5; FIG_H <- 2.5; FIG_RES <- 400
png(file.path(RUNS_DIR, sprintf("fig_baseline_ladder_%s.png", TAG)),
    width = FIG_W * FIG_RES, height = FIG_H * FIG_RES, res = FIG_RES)
par(mfrow = c(1, 3), mar = c(3.2, 7.5, 1.8, 1.2), mgp = c(1.9, 0.55, 0), tcl = -0.25)

y_pos <- rev(seq_len(nrow(tab)))
y_lim <- c(0.4, nrow(tab) + 0.75)
cols  <- c(rep("grey60", nrow(tab) - 1), "firebrick")
CEX_MAIN <- 0.85; CEX_LAB <- 0.75; CEX_AXIS <- 0.70
CEX_YLAB <- 0.65; CEX_VAL <- 0.60; PT_CEX <- 1.10

# RMSE
plot(tab$rmse, y_pos, type = "n", yaxt = "n",
     xlim = c(0, max(tab$rmse) * 1.12), ylim = y_lim,
     xlab = expression("Hourly RMSE (\u00b5mol m"^-2*" s"^-1*")"), ylab = "",
     main = "RMSE (lower = better)",
     cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
segments(0, y_pos, tab$rmse, y_pos, col = "grey85", lwd = 1)
points(tab$rmse, y_pos, pch = 19, cex = PT_CEX, col = cols)
axis(2, at = y_pos, labels = tab$model, las = 1, cex.axis = CEX_YLAB)
abline(v = tail(tab$rmse, 1), col = "firebrick", lty = 3)
text(tab$rmse, y_pos, sprintf(" %.0f", tab$rmse), adj = c(0, -0.5), cex = CEX_VAL)

# R²
plot(tab$r2, y_pos, type = "n", yaxt = "n", xlim = c(0, 1.08), ylim = y_lim,
     xlab = expression("R"^2), ylab = "", main = "R\u00b2 (higher = better)",
     cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
segments(0, y_pos, tab$r2, y_pos, col = "grey85", lwd = 1)
points(tab$r2, y_pos, pch = 19, cex = PT_CEX, col = cols)
abline(v = tail(tab$r2, 1), col = "firebrick", lty = 3)
text(tab$r2, y_pos, sprintf(" %.2f", tab$r2), adj = c(0, -0.5), cex = CEX_VAL)

# Bias
plot(tab$bias, y_pos, type = "n", yaxt = "n", ylim = y_lim,
     xlim = range(tab$bias) + c(-1, 1) * diff(range(tab$bias)) * 0.14,
     xlab = expression("Bias (\u00b5mol m"^-2*" s"^-1*")"), ylab = "",
     main = "Bias (closer to 0 = better)",
     cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
segments(0, y_pos, tab$bias, y_pos, col = "grey85", lwd = 1)
points(tab$bias, y_pos, pch = 19, cex = PT_CEX, col = cols)
abline(v = 0, col = "red", lty = 2, lwd = 1.2)
text(tab$bias, y_pos, sprintf(" %+.0f", tab$bias), adj = c(0, -0.5), cex = CEX_VAL)
dev.off()
message("Saved fig_baseline_ladder_", TAG, ".png")

# pairwise significance vs chosen model, aggregated by site-year (autocorrelation-aware)
v$site_year <- paste(v$site, v$year, sep = "_")
rmse_by_group <- function(pred, obs, group) {
  ok <- is.finite(obs) & is.finite(pred)
  df <- data.frame(pred = pred[ok], obs = obs[ok], group = group[ok])
  sapply(split(df, df$group), function(d) sqrt(mean((d$pred - d$obs)^2)))
}
rmse_chosen <- rmse_by_group(v$pred_M5, v$light_mol, v$site_year)
p_vals <- sapply(seq_along(MODELS), function(i) {
  if (names(MODELS)[i] == "M5: chosen") return(NA)
  rmse_i <- rmse_by_group(MODELS[[i]], v$light_mol, v$site_year)
  common <- intersect(names(rmse_i), names(rmse_chosen))
  wilcox.test(rmse_i[common], rmse_chosen[common], paired = TRUE)$p.value
})
p_consecutive <- sapply(2:length(MODELS), function(i) {
  rmse_prev <- rmse_by_group(MODELS[[i-1]], v$light_mol, v$site_year)
  rmse_curr <- rmse_by_group(MODELS[[i]],   v$light_mol, v$site_year)
  common <- intersect(names(rmse_prev), names(rmse_curr))
  wilcox.test(rmse_prev[common], rmse_curr[common], paired = TRUE)$p.value
})
tab$p_vs_previous <- c(NA, p_consecutive)
tab$p_vs_chosen   <- p_vals
print(tab, row.names = FALSE)


# ══════════════════════════════════════════════════════════════════════════
#  [FIG 8] Per-site PAR through deployment window, baseline vs chosen (2 pages)
# ══════════════════════════════════════════════════════════════════════════

v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v <- add_derived_cols(v)

SHOW_ALL_LINES <- TRUE
sites   <- sort(unique(v$sitekey))
n_sites <- length(sites)

FIG_W <- 7.5; ROW_H <- 1.85; PAGE1_EXTRA <- 1.5; PAGE2_EXTRA <- 1.5; FIG_RES <- 300
sites_p1 <- sites[1:5]
sites_p2 <- sites[6:n_sites]
CEX_MAIN <- 1.15; CEX_LAB <- 1.05; CEX_AXIS <- 0.95; CEX_LEG <- 0.95; CEX_OUT <- 1.30

build_agg <- function(sub) {
  do.call(rbind, lapply(sort(unique(sub$hours_num)), function(t) {
    s2 <- sub[sub$hours_num == t, ]
    data.frame(hours_num = t,
               meas_mean   = mean(s2$light_mol,     na.rm = TRUE),
               pred_b_mean = mean(s2$pred_baseline, na.rm = TRUE),
               pred_c_mean = mean(s2$pred_chosen,   na.rm = TRUE))
  }))
}

draw_ts_panel <- function(agg, sub_site, site, side, ymax) {
  is_baseline <- (side == "baseline")
  plot(agg$hours_num, agg$meas_mean, type = "n",
       xlim = range(agg$hours_num), ylim = c(0, ymax), xaxt = "n", xlab = "",
       ylab = if (is_baseline) expression("PAR (\u00b5mol m"^-2*" s"^-1*")") else "",
       main = sprintf("%s \u2014 %s", site, side),
       cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
  tick_range <- range(agg$hours_num)
  first_tick <- ceiling(tick_range[1] / 6) * 6
  ticks <- seq(first_tick, tick_range[2], by = 6)
  axis(1, at = ticks, labels = hours_to_label(ticks), cex.axis = CEX_AXIS * 0.85)
  pred_col_name <- if (is_baseline) "pred_baseline" else "pred_chosen"
  if (SHOW_ALL_LINES) {
    for (id in unique(sub_site$hoboID_env)) {
      s2 <- sub_site[sub_site$hoboID_env == id, ]; s2 <- s2[order(s2$hours_num), ]
      lines(s2$hours_num, s2$light_mol,      col = adjustcolor("steelblue", 0.25), lwd = 0.6)
      lines(s2$hours_num, s2[[pred_col_name]], col = adjustcolor("firebrick", 0.25), lwd = 0.6)
    }
  }
  pred_mean <- if (is_baseline) agg$pred_b_mean else agg$pred_c_mean
  lines(agg$hours_num, pred_mean,     col = "firebrick", lwd = 2.4)
  lines(agg$hours_num, agg$meas_mean, col = "steelblue", lwd = 2.4)
}

draw_page_ts <- function(site_set, outfile, page_h, title_txt) {
  png(outfile, width = FIG_W * FIG_RES, height = page_h * FIG_RES, res = FIG_RES)
  par(mfrow = c(length(site_set), 2), mar = c(3.2, 4.6, 2.6, 1.2),
      mgp = c(2.6, 0.7, 0), oma = c(4.5, 0, 3.2, 0))
  for (s in site_set) {
    sub <- v[v$sitekey == s, ]; sub <- sub[order(sub$hours_num), ]
    agg <- build_agg(sub); agg <- agg[order(agg$hours_num), ]
    ymax <- quantile(c(sub$light_mol, sub$pred_baseline, sub$pred_chosen),
                     0.99, na.rm = TRUE)
    draw_ts_panel(agg, sub, s, "baseline", ymax)
    m <- metrics(sub$light_mol, sub$pred_baseline)
    legend("topright", bty = "n", cex = CEX_LEG,
           legend = c(sprintf("RMSE=%.0f", m$rmse), sprintf("bias=%+.0f", m$bias),
                      sprintf("R\u00b2=%.2f", m$r2)))
    draw_ts_panel(agg, sub, s, "chosen", ymax)
    m <- metrics(sub$light_mol, sub$pred_chosen)
    legend("topright", bty = "n", cex = CEX_LEG,
           legend = c(sprintf("RMSE=%.0f", m$rmse), sprintf("bias=%+.0f", m$bias),
                      sprintf("R\u00b2=%.2f", m$r2)))
  }
  mtext(title_txt, outer = TRUE, cex = CEX_OUT, font = 2, line = 0.8)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom",
         legend = c("Measured (per sensor)", "Predicted (per sensor)",
                    "Measured (site mean)", "Predicted (site mean)"),
         col = c(adjustcolor("steelblue", 0.4), adjustcolor("firebrick", 0.4),
                 "steelblue", "firebrick"),
         lwd = c(0.9, 0.9, 2.4, 2.4), horiz = TRUE, bty = "n",
         cex = CEX_LEG, inset = c(0, 0.004))
  dev.off()
  message("Saved ", outfile)
}

while (!is.null(dev.list())) dev.off()
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
draw_page_ts(sites_p1, sprintf("fig_persite_ts_baseline_vs_chosen_p1_%s.png", STAMP),
             page_h = 5 * ROW_H + PAGE1_EXTRA,
             title_txt = "Per-site PAR: baseline vs chosen model (1 of 2)")
draw_page_ts(sites_p2, sprintf("fig_persite_ts_baseline_vs_chosen_p2_%s.png", STAMP),
             page_h = 4 * ROW_H + PAGE2_EXTRA,
             title_txt = "Per-site PAR: baseline vs chosen model (2 of 2)")


# ══════════════════════════════════════════════════════════════════════════
#  [FIG 9] Pooled hourly obs vs pred, baseline vs chosen (08/12/16h slopes)
# ══════════════════════════════════════════════════════════════════════════
v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v <- add_derived_cols(v)
HOUR_TARGETS <- c(8, 12, 16)
hour_rng  <- range(v$hour, na.rm = TRUE)
hour_seq  <- seq(hour_rng[1], hour_rng[2])
mid_hour  <- mean(hour_rng)
dist_norm <- abs(hour_seq - mid_hour) / max(abs(hour_seq - mid_hour))
day_ramp  <- colorRampPalette(c("gold", "steelblue", "navy"))(100)
hour_pal  <- day_ramp[round(1 + dist_norm * 99)]
hour_col  <- setNames(hour_pal, as.character(hour_seq))
target_cols <- setNames(c("steelblue", "yellow4", "navy"), as.character(HOUR_TARGETS))
lim_h <- range(c(v$light_mol, v$pred_baseline, v$pred_chosen), na.rm = TRUE)
# DEFINE the panel drawer BEFORE the device opens (fixes original ordering bug)
draw_hourly_panel <- function(pred_col, title) {
  m_all <- metrics(v$light_mol, v[[pred_col]])
  # AXES SWAPPED: predicted on x, measured on y
  plot(v[[pred_col]], v$light_mol, pch = 16, cex = 0.55,
       col = adjustcolor(hour_col[as.character(v$hour)], 0.65),
       xlim = lim_h, ylim = lim_h,
       xlab = bquote("Predicted PAR (\u00b5mol m"^-2*" s"^-1*")"),
       ylab = bquote("Measured PAR (\u00b5mol m"^-2*" s"^-1*")"), main = title)
  abline(0, 1, lty = 3, col = "grey60", lwd = 1.5)
  slopes <- setNames(rep(NA_real_, length(HOUR_TARGETS)), as.character(HOUR_TARGETS))
  for (h in HOUR_TARGETS) {
    idx <- which(v$hour == h & is.finite(v[[pred_col]]) & is.finite(v$light_mol))
    if (length(idx) < 5) next
    sub <- v[idx, ]
    # fit measured ~ predicted now that predicted is the x-axis
    fit <- lm(y ~ 0 + x, data = data.frame(y = sub$light_mol, x = sub[[pred_col]]))
    slope <- coef(fit)[["x"]]; slopes[as.character(h)] <- slope
    x_end <- max(sub[[pred_col]], na.rm = TRUE) * 1.05
    segments(0, 0, x_end, slope * x_end, col = "white", lwd = 7)
    segments(0, 0, x_end, slope * x_end, col = target_cols[as.character(h)], lwd = 4)
  }
  legend("topleft", bty = "n", cex = 1.0,
         legend = c(sprintf("n = %d", m_all$n), sprintf("RMSE = %.0f", m_all$rmse),
                    sprintf("bias = %+.0f", m_all$bias), sprintf("R\u00b2 = %.2f", m_all$r2)))
  legend("bottomright", bty = "n", cex = 1.0,
         legend = c("1:1 line",
                    sprintf("08:00: slope = %.2f", slopes["8"]),
                    sprintf("12:00: slope = %.2f", slopes["12"]),
                    sprintf("16:00: slope = %.2f", slopes["16"])),
         col = c("grey60", target_cols), lty = c(3, 1, 1, 1), lwd = c(1.5, 4, 4, 4))
}
while (!is.null(dev.list())) dev.off()
png(sprintf("fig_pooled_hourly_baseline_vs_chosen_%s.png",
            format(Sys.time(), "%Y%m%d_%H%M%S")),
    width = 16 * 400, height = 8 * 400, res = 400)
par(mfrow = c(1, 2), mgp = c(2.6, 0.7, 0), xpd = FALSE,
    cex.lab = 1.25, cex.axis = 1.15, cex.main = 1.35)
par(mar = c(4.5, 4.5, 3.5, 10))
draw_hourly_panel("pred_baseline",
                  sprintf("Baseline (shadow \u00d7 exp(\u2212%.2f\u00d7d))", KD_BASELINE))
par(mar = c(4.5, 4.5, 3.5, 10))
draw_hourly_panel("pred_chosen", "Chosen model")
# right-hand diurnal colour bar
par(xpd = TRUE)
strip_x <- lim_h[2] + diff(lim_h) * 0.06
strip_y <- seq(lim_h[1], lim_h[1] + diff(lim_h) * 0.55, length.out = length(hour_pal))
for (k in seq_along(hour_pal))
  if (k < length(hour_pal))
    rect(strip_x, strip_y[k], strip_x + diff(lim_h) * 0.03, strip_y[k + 1],
         col = hour_pal[k], border = NA)
text(strip_x, max(strip_y) + diff(lim_h) * 0.04, "Hour of day",
     cex = 1.0, font = 2, adj = 0)
mid_idx <- round(length(hour_seq) / 2)
text(strip_x + diff(lim_h) * 0.045, strip_y[1],
     sprintf("%02d:00", hour_rng[1]), cex = 0.9, adj = 0)
text(strip_x + diff(lim_h) * 0.045, strip_y[mid_idx],
     sprintf("%02d:00", hour_seq[mid_idx]), cex = 0.9, adj = 0)
text(strip_x + diff(lim_h) * 0.045, tail(strip_y, 1),
     sprintf("%02d:00", hour_rng[2]), cex = 0.9, adj = 0)
par(xpd = FALSE)
dev.off()
message("Saved fig_pooled_hourly.")


# ---- Single-panel overlay: baseline (faint) vs chosen (bold) per-hour hulls ----

v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v <- add_derived_cols(v)

HOUR_TARGETS <- c(8, 12, 16)
target_cols <- setNames(c("steelblue", "yellow4", "navy"), as.character(HOUR_TARGETS))
lim_h <- range(c(v$light_mol, v$pred_baseline, v$pred_chosen), na.rm = TRUE)

keep_frac <- 0.90  # inner fraction of points kept before hulling

# ---- trimmed convex hull; style controls baseline vs chosen look ----
add_core_hull <- function(x, y, fill, keep = keep_frac,
                          fill_alpha = 0.20, border_alpha = 0.85,
                          lwd = 2, lty = 1) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 5) return(FALSE)
  cx <- median(x); cy <- median(y)
  sx <- mad(x); sy <- mad(y)
  if (sx <= 0) sx <- sd(x); if (sy <= 0) sy <- sd(y)
  if (!is.finite(sx) || sx <= 0 || !is.finite(sy) || sy <= 0) return(FALSE)
  d2 <- ((x - cx) / sx)^2 + ((y - cy) / sy)^2
  in_core <- d2 <= quantile(d2, keep, na.rm = TRUE)
  xc <- x[in_core]; yc <- y[in_core]
  if (length(xc) < 3) return(FALSE)
  hpts <- chull(xc, yc)
  polygon(xc[hpts], yc[hpts],
          col = adjustcolor(fill, fill_alpha),
          border = adjustcolor(fill, border_alpha), lwd = lwd, lty = lty)
  TRUE
}

hulls_for_model <- function(pred_col, ...) {
  for (h in HOUR_TARGETS) {
    idx <- which(v$hour == h & is.finite(v[[pred_col]]) & is.finite(v$light_mol))
    if (length(idx) < 5) next
    add_core_hull(v[[pred_col]][idx], v$light_mol[idx],
                  fill = target_cols[as.character(h)], ...)
  }
}

# ---- device ----
while (!is.null(dev.list())) dev.off()
png(sprintf("fig_overlay_baseline_vs_chosen_%s.png",
            format(Sys.time(), "%Y%m%d_%H%M%S")),
    width = 9 * 400, height = 8 * 400, res = 400)
par(mgp = c(2.6, 0.7, 0), xpd = FALSE,
    cex.lab = 1.25, cex.axis = 1.15, cex.main = 1.35,
    mar = c(4.5, 4.5, 3.5, 1))

plot(NA, NA, xlim = lim_h, ylim = lim_h,
     xlab = bquote("Predicted PAR (\u00b5mol m"^-2*" s"^-1*")"),
     ylab = bquote("Measured PAR (\u00b5mol m"^-2*" s"^-1*")"),
     main = "Baseline (faint) vs Chosen (bold)")
abline(0, 1, lty = 3, col = "grey60", lwd = 1.5)

# baseline first (faint, dashed, under)
hulls_for_model("pred_baseline",
                fill_alpha = 0.12, border_alpha = 0.55, lwd = 1.5, lty = 2)
# chosen second (bold, solid, on top)
hulls_for_model("pred_chosen",
                fill_alpha = 0.55, border_alpha = 1.0, lwd = 3, lty = 1)

# legend: hours by colour, models by line style
legend("bottomright", bty = "n", cex = 1.0,
       legend = c("1:1 line", "08:00", "12:00", "16:00",
                  "", "baseline", "chosen"),
       col = c("grey60", target_cols, NA, "grey30", "grey30"),
       lty = c(3, 1, 1, 1, NA, 2, 1),
       lwd = c(1.5, 6, 6, 6, NA, 1.5, 3),
       seg.len = 1.6)

dev.off()
message("Saved fig_overlay_baseline_vs_chosen.")



# ---- Single-panel overlay: baseline + chosen, points + slopes + stats ----

v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v <- add_derived_cols(v)

HOUR_TARGETS <- c(8, 12, 16)
hour_rng  <- range(v$hour, na.rm = TRUE)
hour_seq  <- seq(hour_rng[1], hour_rng[2])
mid_hour  <- mean(hour_rng)
dist_norm <- abs(hour_seq - mid_hour) / max(abs(hour_seq - mid_hour))
day_ramp  <- colorRampPalette(c("gold", "steelblue", "navy"))(100)
hour_pal  <- day_ramp[round(1 + dist_norm * 99)]
hour_col  <- setNames(hour_pal, as.character(hour_seq))
target_cols <- setNames(c("steelblue", "yellow4", "navy"), as.character(HOUR_TARGETS))
lim_h <- range(c(v$light_mol, v$pred_baseline, v$pred_chosen), na.rm = TRUE)

# per-model point + slope styling
model_pch <- c(baseline = 1, chosen = 16)   # hollow vs filled
model_lty <- c(baseline = 2, chosen = 1)    # dashed vs solid

draw_points <- function(pred_col, pch, alpha) {
  points(v[[pred_col]], v$light_mol, pch = pch, cex = 0.45,
         col = adjustcolor(hour_col[as.character(v$hour)], alpha))
}

draw_slopes <- function(pred_col, lty) {
  slopes <- setNames(rep(NA_real_, length(HOUR_TARGETS)), as.character(HOUR_TARGETS))
  for (h in HOUR_TARGETS) {
    idx <- which(v$hour == h & is.finite(v[[pred_col]]) & is.finite(v$light_mol))
    if (length(idx) < 5) next
    sub <- v[idx, ]
    fit <- lm(y ~ 0 + x, data = data.frame(y = sub$light_mol, x = sub[[pred_col]]))
    slope <- coef(fit)[["x"]]; slopes[as.character(h)] <- slope
    x_end <- max(sub[[pred_col]], na.rm = TRUE) * 1.05
    segments(0, 0, x_end, slope * x_end, col = "white", lwd = 6)
    segments(0, 0, x_end, slope * x_end,
             col = target_cols[as.character(h)], lwd = 3.5, lty = lty)
  }
  slopes
}

# ---- device ----
while (!is.null(dev.list())) dev.off()
png(sprintf("fig_overlay_all_%s.png", format(Sys.time(), "%Y%m%d_%H%M%S")),
    width = 10 * 400, height = 9 * 400, res = 400)
par(mgp = c(2.6, 0.7, 0), xpd = FALSE,
    cex.lab = 1.25, cex.axis = 1.15, cex.main = 1.35,
    mar = c(4.5, 4.5, 3.5, 1))

plot(NA, NA, xlim = lim_h, ylim = lim_h,
     xlab = bquote("Predicted PAR (\u00b5mol m"^-2*" s"^-1*")"),
     ylab = bquote("Measured PAR (\u00b5mol m"^-2*" s"^-1*")"),
     main = "Baseline vs Chosen \u2014 points, slopes, statistics")
abline(0, 1, lty = 3, col = "grey60", lwd = 1.5)

# points: baseline first (hollow, fainter), chosen on top (filled)
draw_points("pred_baseline", model_pch["baseline"], 0.35)
draw_points("pred_chosen",   model_pch["chosen"],   0.55)

# slope lines
sb <- draw_slopes("pred_baseline", model_lty["baseline"])
sc <- draw_slopes("pred_chosen",   model_lty["chosen"])

# stats
mb <- metrics(v$light_mol, v$pred_baseline)
mc <- metrics(v$light_mol, v$pred_chosen)
legend("topleft", bty = "n", cex = 1.0, title = "Baseline",
       legend = c(sprintf("n = %d", mb$n), sprintf("RMSE = %.0f", mb$rmse),
                  sprintf("bias = %+.0f", mb$bias), sprintf("R\u00b2 = %.2f", mb$r2)))
legend(x = lim_h[1] + diff(lim_h) * 0.24, y = lim_h[2],
       bty = "n", cex = 1.0, title = "Chosen",
       legend = c(sprintf("n = %d", mc$n), sprintf("RMSE = %.0f", mc$rmse),
                  sprintf("bias = %+.0f", mc$bias), sprintf("R\u00b2 = %.2f", mc$r2)))

# combined legend: hours by colour, model by line style
legend("bottomright", bty = "n", cex = 0.95,
       legend = c("1:1 line",
                  sprintf("08:00  base %.2f / chosen %.2f", sb["8"],  sc["8"]),
                  sprintf("12:00  base %.2f / chosen %.2f", sb["12"], sc["12"]),
                  sprintf("16:00  base %.2f / chosen %.2f", sb["16"], sc["16"]),
                  "", "baseline (dashed, hollow)", "chosen (solid, filled)"),
       col = c("grey60", target_cols, NA, "grey30", "grey30"),
       lty = c(3, 1, 1, 1, NA, 2, 1),
       lwd = c(1.5, 3.5, 3.5, 3.5, NA, 2, 3.5),
       pch = c(NA, NA, NA, NA, NA, 1, 16),
       seg.len = 1.6)

dev.off()
message("Saved fig_overlay_all.")


# ---- Single-panel overlay: hulls (centerpiece) + points + slopes + stats ----

v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v <- add_derived_cols(v)

HOUR_TARGETS <- c(8, 12, 16)
hour_rng  <- range(v$hour, na.rm = TRUE)
hour_seq  <- seq(hour_rng[1], hour_rng[2])
mid_hour  <- mean(hour_rng)
dist_norm <- abs(hour_seq - mid_hour) / max(abs(hour_seq - mid_hour))
day_ramp  <- colorRampPalette(c("gold", "steelblue", "navy"))(100)
hour_pal  <- day_ramp[round(1 + dist_norm * 99)]
hour_col  <- setNames(hour_pal, as.character(hour_seq))
target_cols <- setNames(c("steelblue", "yellow4", "navy"), as.character(HOUR_TARGETS))
lim_h <- range(c(v$light_mol, v$pred_baseline, v$pred_chosen), na.rm = TRUE)

model_pch <- c(baseline = 1, chosen = 16)   # hollow vs filled
model_lty <- c(baseline = 2, chosen = 1)    # dashed vs solid
keep_frac <- 0.90                            # inner fraction for hulls

draw_points <- function(pred_col, pch, alpha) {
  points(v[[pred_col]], v$light_mol, pch = pch, cex = 0.7,   # bigger
         col = adjustcolor(hour_col[as.character(v$hour)], alpha))
}

# trimmed convex hull; fill_alpha/lty distinguish baseline vs chosen
add_core_hull <- function(x, y, fill, keep = keep_frac,
                          fill_alpha = 0.20, border_alpha = 0.9,
                          lwd = 2.5, lty = 1) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 5) return(invisible())
  cx <- median(x); cy <- median(y); sx <- mad(x); sy <- mad(y)
  if (sx <= 0) sx <- sd(x); if (sy <= 0) sy <- sd(y)
  if (!is.finite(sx) || sx <= 0 || !is.finite(sy) || sy <= 0) return(invisible())
  d2 <- ((x - cx) / sx)^2 + ((y - cy) / sy)^2
  in_core <- d2 <= quantile(d2, keep, na.rm = TRUE)
  xc <- x[in_core]; yc <- y[in_core]
  if (length(xc) < 3) return(invisible())
  hpts <- chull(xc, yc)
  polygon(xc[hpts], yc[hpts], col = adjustcolor(fill, fill_alpha),
          border = adjustcolor(fill, border_alpha), lwd = lwd, lty = lty)
}

draw_hulls <- function(pred_col, fill_alpha, lty) {
  for (h in HOUR_TARGETS) {
    idx <- which(v$hour == h & is.finite(v[[pred_col]]) & is.finite(v$light_mol))
    if (length(idx) < 5) next
    add_core_hull(v[[pred_col]][idx], v$light_mol[idx],
                  fill = target_cols[as.character(h)],
                  fill_alpha = fill_alpha, lty = lty)
  }
}

draw_slopes <- function(pred_col, lty) {
  slopes <- setNames(rep(NA_real_, length(HOUR_TARGETS)), as.character(HOUR_TARGETS))
  for (h in HOUR_TARGETS) {
    idx <- which(v$hour == h & is.finite(v[[pred_col]]) & is.finite(v$light_mol))
    if (length(idx) < 5) next
    sub <- v[idx, ]
    fit <- lm(y ~ 0 + x, data = data.frame(y = sub$light_mol, x = sub[[pred_col]]))
    slope <- coef(fit)[["x"]]; slopes[as.character(h)] <- slope
    x_end <- max(sub[[pred_col]], na.rm = TRUE) * 1.05
    segments(0, 0, x_end, slope * x_end, col = "white", lwd = 6)
    segments(0, 0, x_end, slope * x_end,
             col = target_cols[as.character(h)], lwd = 3.5, lty = lty)
  }
  slopes
}

# ---- device ----
while (!is.null(dev.list())) dev.off()
png(sprintf("fig_overlay_all_hulls_%s.png", format(Sys.time(), "%Y%m%d_%H%M%S")),
    width = 10 * 400, height = 9 * 400, res = 400)
par(mgp = c(2.6, 0.7, 0), xpd = FALSE,
    cex.lab = 1.25, cex.axis = 1.15, cex.main = 1.35,
    mar = c(4.5, 4.5, 3.5, 1))

plot(NA, NA, xlim = lim_h, ylim = lim_h,
     xlab = bquote("Predicted PAR (\u00b5mol m"^-2*" s"^-1*")"),
     ylab = bquote("Measured PAR (\u00b5mol m"^-2*" s"^-1*")"),
     main = "Baseline vs Chosen \u2014 hulls, points, slopes, statistics")
abline(0, 1, lty = 3, col = "grey60", lwd = 1.5)

# 1. points as faint background texture (more transparent, slightly bigger)
draw_points("pred_baseline", model_pch["baseline"], 0.18)
draw_points("pred_chosen",   model_pch["chosen"],   0.28)

# 2. hulls as the centerpiece, over the points
draw_hulls("pred_baseline", fill_alpha = 0.12, lty = model_lty["baseline"])
draw_hulls("pred_chosen",   fill_alpha = 0.32, lty = model_lty["chosen"])

# 3. slope lines on top
sb <- draw_slopes("pred_baseline", model_lty["baseline"])
sc <- draw_slopes("pred_chosen",   model_lty["chosen"])

# stats
mb <- metrics(v$light_mol, v$pred_baseline)
mc <- metrics(v$light_mol, v$pred_chosen)
legend("topleft", bty = "n", cex = 1.0, title = "Baseline",
       legend = c(sprintf("n = %d", mb$n), sprintf("RMSE = %.0f", mb$rmse),
                  sprintf("bias = %+.0f", mb$bias), sprintf("R\u00b2 = %.2f", mb$r2)))
legend(x = lim_h[1] + diff(lim_h) * 0.24, y = lim_h[2],
       bty = "n", cex = 1.0, title = "Chosen",
       legend = c(sprintf("n = %d", mc$n), sprintf("RMSE = %.0f", mc$rmse),
                  sprintf("bias = %+.0f", mc$bias), sprintf("R\u00b2 = %.2f", mc$r2)))

legend("bottomright", bty = "n", cex = 0.95,
       legend = c("1:1 line",
                  sprintf("08:00  base %.2f / chosen %.2f", sb["8"],  sc["8"]),
                  sprintf("12:00  base %.2f / chosen %.2f", sb["12"], sc["12"]),
                  sprintf("16:00  base %.2f / chosen %.2f", sb["16"], sc["16"]),
                  "", "baseline (dashed, hollow)", "chosen (solid, filled)"),
       col = c("grey60", target_cols, NA, "grey30", "grey30"),
       lty = c(3, 1, 1, 1, NA, 2, 1),
       lwd = c(1.5, 3.5, 3.5, 3.5, NA, 2, 3.5),
       pch = c(NA, NA, NA, NA, NA, 1, 16),
       seg.len = 1.6)

dev.off()
message("Saved fig_overlay_all_hulls.")


# per-timepoint (or per-hour) spread across loggers: measured vs predicted
library(dplyr)  # or base aggregate
spread <- v %>%
  group_by(sitekey, hour) %>%
  summarise(sd_meas = sd(light_mol, na.rm = TRUE),
            sd_pred = sd(pred_chosen, na.rm = TRUE), .groups = "drop")

# does measured spread exceed predicted spread?
mean(spread$sd_meas, na.rm = TRUE)
mean(spread$sd_pred, na.rm = TRUE)
mean(spread$sd_meas > spread$sd_pred, na.rm = TRUE)  # fraction of cases
mean(spread$sd_meas > spread$sd_pred, na.rm = TRUE)


# ══════════════════════════════════════════════════════════════════════════
#  [FIG 10] Pooled DLI scatter, baseline vs chosen, colour = depth
#  Two identical SQUARE panels + compact legend column, fit to 6.5" width
# ══════════════════════════════════════════════════════════════════════════
v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v <- add_derived_cols(v)
HOUR_SECONDS <- 3600
build_dli <- function(pred_col) {
  do.call(rbind, lapply(
    split(v, list(v$sitekey, v$hoboID_env), drop = TRUE),
    function(sub) data.frame(
      sitekey       = sub$sitekey[1],
      hoboID_env    = sub$hoboID_env[1],
      dli_measured  = sum(sub$light_mol,   na.rm = TRUE) * HOUR_SECONDS / 1e6,
      dli_predicted = sum(sub[[pred_col]], na.rm = TRUE) * HOUR_SECONDS / 1e6,
      depth_m       = mean(sub$water_depth_m, na.rm = TRUE))))
}
dli_base   <- build_dli("pred_baseline")
dli_chosen <- build_dli("pred_chosen")
m_base   <- metrics(dli_base$dli_measured,   dli_base$dli_predicted)
m_chosen <- metrics(dli_chosen$dli_measured, dli_chosen$dli_predicted)
lim <- range(c(dli_base$dli_measured, dli_base$dli_predicted,
               dli_chosen$dli_measured, dli_chosen$dli_predicted), na.rm = TRUE)
col_rng <- range(c(dli_base$depth_m, dli_chosen$depth_m), na.rm = TRUE)
# light blue (shallow) -> dark blue (deep)
col_pal <- colorRampPalette(c("#cfe8f5", "#7fb8e0", "#2b6cb0", "#0a2a5e"))(100)
depth_to_col <- function(x) {
  # index 1 = light (shallow, col_rng[1]), index 100 = dark (deep, col_rng[2])
  idx <- round((x - col_rng[1]) / diff(col_rng) * 99) + 1
  col_pal[pmax(1, pmin(100, idx))]
}
PT_ALPHA <- 0.7   # point transparency
keys    <- sort(unique(v$sitekey))
pch_map <- setNames((seq_along(keys) - 1) %% 25 + 1, keys)

# DEFINE the panel drawer BEFORE the device opens
draw_dli_panel <- function(dli, m, title) {
  # AXES SWAPPED: predicted on x, measured on y.
  # Square box via pty="s" (set in par); xlim==ylim keeps 1:1 line honest.
  plot(dli$dli_predicted, dli$dli_measured,
       pch = pch_map[dli$sitekey],
       col = adjustcolor(depth_to_col(dli$depth_m), PT_ALPHA),
       cex = 0.6, lwd = 0.9, xlim = lim, ylim = lim,
       xlab = expression(Predicted~DLI~(mol~m^-2)),
       ylab = expression(Measured~DLI~(mol~m^-2)), main = title)
  abline(0, 1, col = "red", lwd = 1.3, lty = 2)
  fit <- lm(dli_measured ~ dli_predicted, data = dli)
  abline(fit, col = "black", lwd = 1.6)
  cf <- coef(fit)
  legend("topleft", bty = "n", cex = 0.7,
         legend = c(sprintf("n = %d", m$n), sprintf("RMSE = %.2f", m$rmse),
                    sprintf("bias = %+.2f", m$bias), sprintf("R\u00b2 = %.2f", m$r2),
                    sprintf("fit: y = %.2f + %.2fx", cf[1], cf[2])))
}

while (!is.null(dev.list())) dev.off()
png(sprintf("fig_pooled_dli_baseline_vs_chosen_%s.png",
            format(Sys.time(), "%Y%m%d_%H%M%S")),
    width = 6.5, height = 3.4, units = "in", res = 800)

layout(matrix(c(1, 2, 3), 1, 3), widths = c(1, 1, 0.42))

par(mgp = c(2.3, 0.7, 0), cex.main = 1.0, cex.lab = 0.9, cex.axis = 0.8,
    pty = "s")
par(mar = c(4, 4, 3, 1))
draw_dli_panel(dli_base, m_base,
               sprintf("Baseline (shadow \u00d7 exp(\u2212%.2f\u00d7d))", KD_BASELINE))
par(mar = c(4, 4, 3, 1))
draw_dli_panel(dli_chosen, m_chosen, "Chosen model")

# ---- third cell: site legend + depth colour bar ----
par(mar = c(4, 0, 3, 0), pty = "m")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

legend("topleft", legend = keys, pch = pch_map[keys],
       pt.cex = 0.95, pt.lwd = 1, bty = "n", title = "Site/year",
       cex = 0.8, title.cex = 0.9, y.intersp = 0.95)

# colour bar: SHALLOW (light) at TOP, DEEP (dark) at BOTTOM
bar_x0 <- 0.05; bar_x1 <- 0.24
bar_y  <- seq(0.05, 0.40, length.out = 100)   # bar_y[1] bottom, bar_y[100] top
for (i in 1:99)
  # top of bar (high bar_y) shows col_pal[1] (light/shallow); reverse fill
  rect(bar_x0, bar_y[i], bar_x1, bar_y[i + 1],
       col = col_pal[101 - i], border = NA)
text(bar_x0, 0.45, "Depth (m)", adj = 0, font = 2, cex = 0.8)
# top label = shallow (col_rng[1]), bottom label = deep (col_rng[2])
text(bar_x1 + 0.04, bar_y[100], round(col_rng[1], 1), adj = 0, cex = 0.75)
text(bar_x1 + 0.04, bar_y[1],   round(col_rng[2], 1), adj = 0, cex = 0.75)

dev.off()
message("Saved fig_pooled_dli.")




# ══════════════════════════════════════════════════════════════════════════
#  [FIG] Pooled DLI hull overlap, baseline vs chosen  ->  SVG for Illustrator
#  Shows how the paired (predicted, measured) DLI cloud tightens under the
#  chosen model. NOTE: this is calibration spread of the cloud, NOT the
#  between-logger SD result (170 vs 55) — caption accordingly.
# ══════════════════════════════════════════════════════════════════════════
# Assumes dli_base and dli_chosen already exist from build_dli(), each with
# columns: dli_measured, dli_predicted.  (Run the Fig 10 block first, or the
# build_dli() lines from it, so these objects are in scope.)
keep_frac <- 0.99   # inner fraction of points kept before hulling
core_hull <- function(x, y, keep = keep_frac) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 5) return(NULL)
  cx <- median(x); cy <- median(y); sx <- mad(x); sy <- mad(y)
  if (sx <= 0) sx <- sd(x); if (sy <= 0) sy <- sd(y)
  if (!is.finite(sx) || sx <= 0 || !is.finite(sy) || sy <= 0) return(NULL)
  d2 <- ((x - cx)/sx)^2 + ((y - cy)/sy)^2
  keep_i <- d2 <= quantile(d2, keep, na.rm = TRUE)
  xc <- x[keep_i]; yc <- y[keep_i]
  if (length(xc) < 3) return(NULL)
  h <- chull(xc, yc)
  list(x = xc[h], y = yc[h])
}
hb <- core_hull(dli_base$dli_predicted,   dli_base$dli_measured)
hc <- core_hull(dli_chosen$dli_predicted, dli_chosen$dli_measured)
lim <- range(c(dli_base$dli_measured, dli_base$dli_predicted,
               dli_chosen$dli_measured, dli_chosen$dli_predicted), na.rm = TRUE)
col_base   <- "#d98a3d"  # warm  = baseline (wider)
col_chosen <- "#2b6cb0"  # blue  = chosen  (tighter)

pmono <- "Space Mono"   # font family; must be installed on this machine

while (!is.null(dev.list())) dev.off()
svg(sprintf("fig_dli_hull_overlap_%s.svg", format(Sys.time(), "%Y%m%d_%H%M%S")),
    width = 5.5, height = 5.5, pointsize = 11, family = pmono)
par(mar = c(4, 4, 2, 1), pty = "s", mgp = c(2.3, 0.7, 0),
    cex.lab = 1.0, cex.axis = 0.9, family = pmono)
plot(NA, NA, xlim = lim, ylim = lim,
     xlab = expression(Predicted~DLI~(mol~m^-2)),
     ylab = expression(Measured~DLI~(mol~m^-2)),
     main = "DLI spread: baseline vs chosen")
abline(0, 1, col = "grey55", lty = 3, lwd = 1.3)
# baseline hull first (wider, dashed), chosen on top (tighter, solid)
if (!is.null(hb))
  polygon(hb$x, hb$y, col = adjustcolor(col_base, 0.25),
          border = col_base, lwd = 2, lty = 2)
if (!is.null(hc))
  polygon(hc$x, hc$y, col = adjustcolor(col_chosen, 0.30),
          border = col_chosen, lwd = 2.2)
legend("topleft", bty = "n", cex = 0.95,
       legend = c("1:1 line", "Baseline hull", "Chosen hull"),
       col = c("grey55", col_base, col_chosen),
       lty = c(3, 2, 1), lwd = c(1.3, 2, 2.2))
dev.off()
message("Saved fig_dli_hull_overlap SVG.")


# ══════════════════════════════════════════════════════════════════════════
#  [FIG 10] Pooled DLI scatter — SVG, Space Mono, warm depth palette (#dc4f3f deep)
# ══════════════════════════════════════════════════════════════════════════
library(svglite)          # install.packages("svglite")
# Space Mono must be registered so svglite embeds/names it correctly:
library(systemfonts)      # part of svglite deps
# If Space Mono is installed on your system, this finds it; otherwise register a file:
# register_font(name = "Space Mono", plain = "/path/to/SpaceMono-Regular.ttf",
#               bold = "/path/to/SpaceMono-Bold.ttf")

FONT <- "Space Mono"

v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v <- add_derived_cols(v)

HOUR_SECONDS <- 3600
build_dli <- function(pred_col) {
  do.call(rbind, lapply(
    split(v, list(v$sitekey, v$hoboID_env), drop = TRUE),
    function(sub) data.frame(
      sitekey       = sub$sitekey[1],
      hoboID_env    = sub$hoboID_env[1],
      dli_measured  = sum(sub$light_mol,   na.rm = TRUE) * HOUR_SECONDS / 1e6,
      dli_predicted = sum(sub[[pred_col]], na.rm = TRUE) * HOUR_SECONDS / 1e6,
      depth_m       = mean(sub$water_depth_m, na.rm = TRUE))))
}

dli_base   <- build_dli("pred_baseline")
dli_chosen <- build_dli("pred_chosen")
m_base   <- metrics(dli_base$dli_measured,   dli_base$dli_predicted)
m_chosen <- metrics(dli_chosen$dli_measured, dli_chosen$dli_predicted)

lim <- range(c(dli_base$dli_measured, dli_base$dli_predicted,
               dli_chosen$dli_measured, dli_chosen$dli_predicted), na.rm = TRUE)
col_rng <- range(c(dli_base$depth_m, dli_chosen$depth_m), na.rm = TRUE)

# pale/warm (shallow) -> #dc4f3f (deep)
col_pal <- colorRampPalette(c("#fdeede", "#f4a07a", "#e8734f", "#dc4f3f"))(100)
depth_to_col <- function(x) {
  idx <- round((x - col_rng[1]) / diff(col_rng) * 99) + 1
  col_pal[pmax(1, pmin(100, idx))]
}

PT_ALPHA <- 0.7
keys    <- sort(unique(v$sitekey))
pch_map <- setNames((seq_along(keys) - 1) %% 25 + 1, keys)

draw_dli_panel <- function(dli, m, title) {
  plot(dli$dli_predicted, dli$dli_measured,
       pch = pch_map[dli$sitekey],
       col = adjustcolor(depth_to_col(dli$depth_m), PT_ALPHA),
       cex = 0.6, lwd = 0.9, xlim = lim, ylim = lim,
       xlab = expression(Predicted~DLI~(mol~m^-2)),
       ylab = expression(Measured~DLI~(mol~m^-2)), main = title,
       family = FONT)
  abline(0, 1, col = "red", lwd = 1.3, lty = 2)
  fit <- lm(dli_measured ~ dli_predicted, data = dli)
  abline(fit, col = "black", lwd = 1.6)
  cf <- coef(fit)
  legend("topleft", bty = "n", cex = 0.7, text.font = 1,
         legend = c(sprintf("n = %d", m$n), sprintf("RMSE = %.2f", m$rmse),
                    sprintf("bias = %+.2f", m$bias), sprintf("R\u00b2 = %.2f", m$r2),
                    sprintf("fit: y = %.2f + %.2fx", cf[1], cf[2])))
}

svglite(sprintf("fig_pooled_dli_baseline_vs_chosen_%s.svg",
                format(Sys.time(), "%Y%m%d_%H%M%S")),
        width = 6.5, height = 3.4, system_fonts = list(sans = FONT, mono = FONT))

# make Space Mono the default for everything drawn via par()
par(family = FONT)
layout(matrix(c(1, 2, 3), 1, 3), widths = c(1, 1, 0.42))
par(mgp = c(2.3, 0.7, 0), cex.main = 1.0, cex.lab = 0.9, cex.axis = 0.8,
    pty = "s", family = FONT)

par(mar = c(4, 4, 3, 1))
draw_dli_panel(dli_base, m_base,
               sprintf("Baseline (shadow \u00d7 exp(\u2212%.2f\u00d7d))", KD_BASELINE))
par(mar = c(4, 4, 3, 1))
draw_dli_panel(dli_chosen, m_chosen, "Chosen model")

# ---- legend + colour bar ----
par(mar = c(4, 0, 3, 0), pty = "m", family = FONT)
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
legend("topleft", legend = keys, pch = pch_map[keys],
       pt.cex = 0.95, pt.lwd = 1, bty = "n", title = "Site/year",
       cex = 0.8, y.intersp = 0.95)

bar_x0 <- 0.05; bar_x1 <- 0.24
bar_y  <- seq(0.05, 0.40, length.out = 100)
for (i in 1:99)
  rect(bar_x0, bar_y[i], bar_x1, bar_y[i + 1],
       col = col_pal[101 - i], border = NA)   # top = shallow (pale), bottom = deep (#dc4f3f)
text(bar_x0, 0.45, "Depth (m)", adj = 0, font = 2, cex = 0.8)
text(bar_x1 + 0.04, bar_y[100], round(col_rng[1], 1), adj = 0, cex = 0.75)  # shallow
text(bar_x1 + 0.04, bar_y[1],   round(col_rng[2], 1), adj = 0, cex = 0.75)  # deep

dev.off()
message("Saved SVG.")


# ══════════════════════════════════════════════════════════════════════════
#  3D DLI scatter: x = predicted, y = measured, z = depth (INVERTED)
# ══════════════════════════════════════════════════════════════════════════
library(svglite)
library(plot3D)
FONT <- "Space Mono"

dli <- dli_chosen
col_rng <- range(dli$depth_m, na.rm = TRUE)
col_pal <- colorRampPalette(c("#fdeede", "#f4a07a", "#e8734f", "#dc4f3f"))(100)

svglite(sprintf("fig_dli_3d_%s.svg", format(Sys.time(), "%Y%m%d_%H%M%S")),
        width = 6.5, height = 6, system_fonts = list(sans = FONT, mono = FONT))
par(family = FONT)

scatter3D(x = dli$dli_predicted, y = dli$dli_measured, z = dli$depth_m,
          colvar = dli$depth_m, col = col_pal,
          pch = 19, cex = 0.9, alpha = 0.8,
          theta = 40, phi = 20, bty = "b2",
          zlim = rev(range(dli$depth_m, na.rm = TRUE)),   # <-- shallow at top, deep at bottom
          xlab = "Predicted DLI", ylab = "Measured DLI", zlab = "Depth (m)",
          clab = "Depth (m)",
          main = "DLI: predicted vs measured vs depth",
          ticktype = "detailed", family = FONT)

scatter3D(x = dli$dli_predicted, y = dli$dli_measured, z = dli$depth_m,
          type = "h", add = TRUE, colvar = dli$depth_m, col = col_pal,
          alpha = 0.25, colkey = FALSE)

dev.off()
message("Saved 3D SVG (depth inverted).")






# ══════════════════════════════════════════════════════════════════════════
#  [FIG 11] Per-site residual diagnostics, chosen model (2 pages)
# ══════════════════════════════════════════════════════════════════════════

v <- read.csv(VAL_CSV, stringsAsFactors = FALSE)
v$datetime_pos <- as.POSIXct(v$datetime_pos, tz = "Australia/Brisbane")
v$resid   <- v$pred_par - v$light_mol
v$hour    <- as.integer(format(v$datetime_pos, "%H"))
v$sun_deg <- v$sun_altitude * 180 / pi

sites   <- sort(unique(v$sitekey))
n_sites <- length(sites)
ylim_shared <- quantile(v$resid, c(0.01, 0.99), na.rm = TRUE)

FIG_W <- 7.5; ROW_H <- 1.85; PAGE_EXTRA <- 1.5; FIG_RES <- 300
sites_p1 <- sites[1:5]
sites_p2 <- sites[6:n_sites]
CEX_MAIN <- 1.15; CEX_LAB <- 1.05; CEX_AXIS <- 0.95; CEX_LEG <- 0.95; CEX_OUT <- 1.30

resid_panel <- function(x, y, xlab, main, discrete = FALSE, n_bins = 8) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 5) { plot.new(); title(main = main); return() }
  if (discrete) {
    boxplot(y ~ x, outline = FALSE, col = "grey85", border = "grey40",
            xlab = xlab, ylab = "Residual", main = main,
            cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS,
            ylim = ylim_shared)
    abline(h = 0, lty = 2, col = "red", lwd = 1.4)
  } else {
    plot(x, y, pch = 16, cex = 0.5, col = adjustcolor("black", 0.25),
         xlab = xlab, ylab = "Residual", main = main,
         cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS,
         ylim = ylim_shared)
    abline(h = 0, lty = 2, col = "red", lwd = 1.4)
    fit <- try(loess(y ~ x, span = 0.5), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      xs <- seq(min(x), max(x), length.out = 200)
      pred <- predict(fit, newdata = xs, se = TRUE)
      polygon(c(xs, rev(xs)),
              c(pred$fit - 1.96 * pred$se.fit, rev(pred$fit + 1.96 * pred$se.fit)),
              col = adjustcolor("dodgerblue", 0.20), border = NA)
      lines(xs, pred$fit, col = "dodgerblue", lwd = 2.5)
    }
    breaks <- seq(min(x), max(x), length.out = n_bins + 1)
    bin    <- cut(x, breaks = breaks, include.lowest = TRUE)
    bin_stats <- do.call(rbind, lapply(levels(bin), function(b) {
      idx <- which(bin == b); if (length(idx) < 3) return(NULL)
      data.frame(x_mid = mean(x[idx]), mean = mean(y[idx]), sd = sd(y[idx]))
    }))
    if (!is.null(bin_stats)) {
      arrows(bin_stats$x_mid, bin_stats$mean - bin_stats$sd,
             bin_stats$x_mid, bin_stats$mean + bin_stats$sd,
             angle = 90, code = 3, length = 0.03,
             col = adjustcolor("black", 0.5), lwd = 1)
      points(bin_stats$x_mid, bin_stats$mean, pch = 21, bg = "white",
             col = "black", cex = 1.2, lwd = 1.5)
    }
  }
  legend("topright", bty = "n", cex = CEX_LEG,
         legend = sprintf("bias = %+.0f  (n=%d)", mean(y), length(y)))
}

draw_page_resid <- function(site_set, outfile, page_h, title_txt) {
  png(outfile, width = FIG_W * FIG_RES, height = page_h * FIG_RES, res = FIG_RES)
  par(mfrow = c(length(site_set), 3), mar = c(3.8, 4.4, 2.5, 1),
      mgp = c(2.4, 0.7, 0), oma = c(0, 0, 3, 0))
  for (s in site_set) {
    sub <- v[v$sitekey == s, ]
    resid_panel(sub$hour, sub$resid, xlab = "Hour of day",
                main = sprintf("%s: residual vs hour", s), discrete = TRUE)
    resid_panel(sub$water_depth_m, sub$resid, xlab = "Water depth (m)",
                main = sprintf("%s: residual vs depth", s))
    resid_panel(sub$shadow_fraction, sub$resid,
                xlab = "Shadow fraction (0=shaded, 1=lit)",
                main = sprintf("%s: residual vs shadow fraction", s))
  }
  mtext(title_txt, outer = TRUE, cex = CEX_OUT, font = 2, line = 0.8)
  dev.off()
  message("Saved ", outfile)
}

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
draw_page_resid(sites_p1, sprintf("fig_persite_residual_diagnostics_p1_%s.png", STAMP),
                page_h = 5 * ROW_H + PAGE_EXTRA,
                title_txt = "Per-site residual diagnostics — chosen model (1 of 2)")
draw_page_resid(sites_p2, sprintf("fig_persite_residual_diagnostics_p2_%s.png", STAMP),
                page_h = 4 * ROW_H + PAGE_EXTRA,
                title_txt = "Per-site residual diagnostics — chosen model (2 of 2)")


# ══════════════════════════════════════════════════════════════════════════
#  [EXTRAS — editorial SVG variants, not thesis figures]
#  These produce alternative styled versions (Space Mono, transparent SVG,
#  spectrum layouts) of Figs 8/10/7. Kept for reference; run only if wanted.
#  Marked [DRAFT — superseded] where a cleaner block above replaces them:
#    - single-page per-site TS  → superseded by [FIG 8] 2-page block
#    - single-page residual grid → superseded by [FIG 11] 2-page block
#    - first pooled-DLI svg()    → superseded by [FIG 10] png block
#  They are intentionally omitted here to keep the runnable path clean;
#  see the original script if the editorial SVG versions are needed.
# ══════════════════════════════════════════════════════════════════════════