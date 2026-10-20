# ══════════════════════════════════════════════════════════════════════════════
#  07_taylor.R: Taylor diagrams of the model ladder and the thesis runs
#
#  Reads : output/predictions.rds (from 03)
#          data/thesis_runs/validation_<tag>.csv, metrics_log.csv (thesis runs)
#  Writes: output/figures/fig_taylor_ladder.png        main text: current models,
#                                                      primary model split into
#                                                      development / new events
#          output/figures/figS5_taylor_with_thesis.png supplement: the same plus
#                                                      every thesis run in grey
#          output/taylor_stats.csv                      the numbers behind both
#
#  [NEW] Claude. Each point is one model scored against the measurements it was
#  validated on: correlation (angle), SD of predictions / SD of measurements
#  (distance from the origin), and centred RMSE / SD of measurements (distance
#  from the "measured" point). Normalising by each set's own measured SD lets
#  runs scored on different events share one diagram, but thesis runs were
#  scored on 5-10 events with errors since fixed (mirrored shading, unstable
#  tide fit), so they are shown for transparency, not as a like-for-like test.
#  A Taylor diagram does not show bias: keep the bias panel of the ladder figure.
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("config.R"))
source(here::here("R", "utils.R"))
suppressPackageStartupMessages({ library(dplyr) })

stage_banner("07  TAYLOR DIAGRAMS")

# events the thesis (exploratory) runs were developed on; everything else in
# SITES was never seen during development
DEV_EVENTS <- c("CB18", "CB19", "TB18", "TB19", "HS19", "N318", "N319", "TM18", "TM19")

taylor_stats <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); o <- obs[ok]; p <- pred[ok]
  so <- stats::sd(o)
  data.frame(n = length(o), r = stats::cor(o, p), sd_ratio = stats::sd(p) / so,
             crmse = sqrt(mean(((p - mean(p)) - (o - mean(o)))^2)) / so,
             bias = mean(p - o))
}

# ── Current models ───────────────────────────────────────────────────────────
preds <- readRDS(need_file(out_path("predictions.rds")))
cur <- preds %>% group_by(model) %>%
  group_modify(~ taylor_stats(.x$light_mol, .x$pred)) %>% ungroup() %>%
  mutate(set = "current", label = model)
split <- preds %>% filter(model == FIG_MAIN_MODEL) %>%
  mutate(set = ifelse(sitekey %in% DEV_EVENTS, "development events", "new events")) %>%
  group_by(set) %>% group_modify(~ taylor_stats(.x$light_mol, .x$pred)) %>% ungroup() %>%
  mutate(model = FIG_MAIN_MODEL, label = paste0(FIG_MAIN_MODEL, ", ", set))

# ── Thesis runs ──────────────────────────────────────────────────────────────
th_dir <- here::here("data", "thesis_runs")
files  <- list.files(th_dir, "^validation_.*\\.csv$", full.names = TRUE)
thesis <- bind_rows(lapply(files, function(f) {
  v <- utils::read.csv(f, stringsAsFactors = FALSE)
  cbind(run_tag = sub("^validation_(.*)\\.csv$", "\\1", basename(f)),
        n_events = length(unique(v$sitekey)),
        taylor_stats(v$light_mol, v$pred_par))
}))
log_f <- file.path(th_dir, "metrics_log.csv")
if (nrow(thesis) > 0 && file.exists(log_f)) {
  lg <- utils::read.csv(log_f, colClasses = c(run_tag = "character"), stringsAsFactors = FALSE)
  thesis <- thesis %>% left_join(lg %>% select(run_tag, kd_method), by = "run_tag")
}
thesis$set <- "thesis run"; thesis$model <- thesis$run_tag; thesis$label <- ""
message(sprintf("  %d current models, %d thesis runs", nrow(cur), nrow(thesis)))

stats_out <- bind_rows(cur, split, thesis)
write.csv(stats_out, out_path("taylor_stats.csv"), row.names = FALSE)

# ── Drawing ──────────────────────────────────────────────────────────────────
LADDER <- c("M0", "M1", "M2s", "M2r", "M3s", "M3r", "M4s", "M4r")
col_of <- function(m) {
  base <- sub("_.*", "", m)
  pal  <- c(M0 = "#7f7f7f", M1 = "#9467bd", M2s = "#8c564b", M2r = "#1f77b4",
            M3s = "#e377c2", M3r = "#17becf", M4s = "#bcbd22", M4r = "#d62728")
  ifelse(base %in% names(pal), pal[base], "#555555")
}

taylor_plot <- function(pts, file, title) {
  pts  <- pts[is.finite(pts$r) & is.finite(pts$sd_ratio), ]
  smax <- max(1.5, ceiling(max(pts$sd_ratio) * 4) / 4 + 0.1)
  png(file, width = 8, height = 7.5, units = "in", res = 200)
  par(mar = c(4.5, 4.5, 3, 1))
  plot(NA, xlim = c(0, smax), ylim = c(0, smax), asp = 1, axes = FALSE,
       xlab = "SD of predictions / SD of measurements", ylab = "", main = title)
  th <- seq(0, pi / 2, length.out = 200)
  for (s in seq(0.25, smax, 0.25)) lines(s * cos(th), s * sin(th), col = "grey88")
  lines(smax * cos(th), smax * sin(th)); segments(0, 0, smax, 0); segments(0, 0, 0, smax)
  axis(1, at = seq(0, smax, 0.25)); axis(2, at = seq(0, smax, 0.25), las = 1)
  for (r in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)) {
    a <- acos(r)
    segments(0, 0, smax * cos(a), smax * sin(a), col = "grey92")
    text(1.04 * smax * cos(a), 1.04 * smax * sin(a), r, cex = 0.7, col = "grey40")
  }
  text(0.80 * smax, 0.80 * smax, "Correlation", srt = -45, col = "grey40", cex = 0.9)
  ph <- seq(0, pi, length.out = 300)                     # centred RMSE contours
  for (e in seq(0.25, 1.5, 0.25)) {
    x <- 1 + e * cos(ph); y <- e * sin(ph); k <- x >= 0 & x^2 + y^2 <= smax^2
    x[!k] <- NA; y[!k] <- NA
    lines(x, y, col = "#a6d3a6", lty = 2)
  }
  points(1, 0, pch = 17, cex = 1.5); text(1, 0, "measured", pos = 1, cex = 0.8)
  xy <- function(d) cbind(d$sd_ratio * d$r, d$sd_ratio * sin(acos(pmin(pmax(d$r, -1), 1))))
  th_pts <- pts[pts$set == "thesis run", ]
  if (nrow(th_pts)) { p <- xy(th_pts); points(p, pch = 1, col = "grey60", cex = 0.9) }
  cu <- pts[pts$set == "current", ]
  if (nrow(cu)) {
    p <- xy(cu); main <- cu$model %in% LADDER
    points(p, pch = ifelse(main, 19, 4), col = col_of(cu$model), cex = ifelse(main, 1.3, 0.9))
    text(p[main, , drop = FALSE], labels = cu$label[main], pos = 4, cex = 0.75,
         col = col_of(cu$model[main]))
  }
  sp <- pts[pts$set %in% c("development events", "new events"), ]
  if (nrow(sp)) {
    p <- xy(sp); points(p, pch = ifelse(sp$set == "new events", 24, 25),
                        bg = col_of(sp$model), col = "black", cex = 1.3)
  }
  legend("topright", bty = "n", cex = 0.8,
         legend = c("current models (labelled)", "sensitivity runs",
                    paste(FIG_MAIN_MODEL, "on development events"),
                    paste(FIG_MAIN_MODEL, "on new events"),
                    if (nrow(th_pts)) "thesis runs (exploratory)" else NULL,
                    "centred RMSE / SD measured"),
         pch = c(19, 4, 25, 24, if (nrow(th_pts)) 1 else NULL, NA),
         lty = c(NA, NA, NA, NA, if (nrow(th_pts)) NA else NULL, 2),
         col = c("grey30", "grey30", "black", "black", if (nrow(th_pts)) "grey60" else NULL, "#a6d3a6"),
         pt.bg = c(NA, NA, col_of(FIG_MAIN_MODEL), col_of(FIG_MAIN_MODEL),
                   if (nrow(th_pts)) NA else NULL, NA))
  invisible(dev.off())
}

taylor_plot(bind_rows(cur, split), out_path("figures", "fig_taylor_ladder.png"),
            "Model ladder, hourly PAR")
if (nrow(thesis))
  taylor_plot(bind_rows(cur, split, thesis), out_path("figures", "figS5_taylor_with_thesis.png"),
              "Current models and thesis runs (exploratory), hourly PAR")

print(as.data.frame(bind_rows(cur, split) %>% select(model, set, n, r, sd_ratio, crmse, bias)),
      digits = 3, row.names = FALSE)
message("07 done.")
