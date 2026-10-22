# ══════════════════════════════════════════════════════════════════════════════
#  utils.R — small helpers shared by every stage
#
#  PROVENANCE: [NEW] unless marked otherwise. Written by Claude to support the
#  split into stages; nothing here changes the science.
# ══════════════════════════════════════════════════════════════════════════════

# [YOURS, EDITED] — your make_key()
#   - [Claude] always two-digit years, so 18 and 2018 both give "TB18"
# ── helper: one canonical site/year key, used everywhere ──────────────────────
# Keys look like "TB18", "CB19" — match the format in kdpar_lookup,
# capture_tide, and the reference loggers so joins/lookups all line up.
make_key <- function(site, year) {
  y <- as.integer(year)
  paste0(site, sprintf("%02d", y %% 100))
}

# Stop with a clear message if an input file is missing.
need_file <- function(path) {
  if (!file.exists(path))
    stop("Missing input: ", path,
         "\nSee README.md for which files go in data/.", call. = FALSE)
  invisible(path)
}

out_path <- function(...) {
  p <- file.path(OUT_DIR, ...)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  p
}

# [NEW] Read a file: from this run's folder if it is there, otherwise from
# BASE_OUT_DIR (the main output/). Lets a labelled sensitivity run reuse the
# stage 01/02 products without copying them. In the main run both are output/.
in_path <- function(...) {
  p <- file.path(OUT_DIR, ...)
  if (file.exists(p)) return(p)
  b <- file.path(BASE_OUT_DIR, ...)
  if (file.exists(b)) return(b)
  out_path(...)
}

git_hash <- function() {
  h <- tryCatch(suppressWarnings(system("git rev-parse --short HEAD", intern = TRUE,
                       ignore.stderr = TRUE)), error = function(e) character())
  if (length(h) == 0) "not a git repository" else h
}

git_dirty <- function() {
  s <- tryCatch(suppressWarnings(system("git status --porcelain", intern = TRUE,
                       ignore.stderr = TRUE)), error = function(e) character())
  length(s) > 0
}

# [NEW] — modelled on your manifest_<RUN_TAG>.txt
# Writes output/manifest_<stage>.txt: when, which commit, every config value,
# and the R session. Called at the end of each stage.
write_manifest <- function(stage, extra = NULL) {
  path <- out_path(sprintf("manifest_%s.txt", stage))
  cfg_names <- ls(envir = globalenv())
  cfg_names <- cfg_names[grepl("^[A-Z][A-Z0-9_]+$", cfg_names)]
  sink(path)
  on.exit(sink())
  cat("STAGE      : ", stage, "\n", sep = "")
  cat("RUN TIME   : ", format(Sys.time(), tz = "UTC"), " UTC\n", sep = "")
  cat("GIT COMMIT : ", git_hash(),
      if (git_dirty()) "  (uncommitted changes present)" else "", "\n\n", sep = "")
  cat("CONFIG\n")
  for (nm in cfg_names) {
    v <- get(nm, envir = globalenv())
    if (is.function(v)) next
    txt <- if (is.data.frame(v)) paste0("<data.frame ", nrow(v), " rows>")
           else paste(format(v), collapse = ", ")
    if (!is.null(names(v)) && !is.data.frame(v))
      txt <- paste(paste0(names(v), "=", format(v)), collapse = ", ")
    cat(sprintf("  %-28s %s\n", nm, txt))
  }
  if (!is.null(extra)) { cat("\nNOTES\n"); cat(paste0("  ", extra, "\n"), sep = "") }
  cat("\nSESSION\n")
  print(sessionInfo())
  sink()
  on.exit()
  log_run(stage, path, headline = attr(extra, "headline"))
  invisible(path)
}

# [NEW] Claude, 22 Oct 2026 — like your metrics_log.csv
#  Every stage run appends one row to output/run_log.csv (one log for every
#  run label), and a dated copy of the stage manifest goes to output/logs/.
#  headline: optional named vector of key results (03 passes its metrics).
log_run <- function(stage, manifest, headline = NULL) {
  logs <- file.path(BASE_OUT_DIR, "logs")
  dir.create(logs, showWarnings = FALSE, recursive = TRUE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  copy  <- file.path(logs, sprintf("%s_%s_%s.txt", stamp, RUN_LABEL, sub("^manifest_", "", stage)))
  file.copy(manifest, copy, overwrite = TRUE)
  row <- data.frame(
    time            = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    run_id          = RUN_ID,
    stage           = stage,
    run_label       = RUN_LABEL,
    output_folder   = sub(paste0("^", here::here(), "/?"), "", OUT_DIR),
    git_commit      = git_hash(),
    uncommitted     = git_dirty(),
    w_ambient       = if (is.finite(W_AMBIENT_FIXED)) as.character(W_AMBIENT_FIXED) else "CV",
    kd_constant     = KD_CONSTANT,
    buffer_m        = BUFFER_M,
    n_events        = length(SITES),
    test_run        = TEST_RUN,
    headline        = if (is.null(headline)) "" else
                        paste(names(headline), headline, sep = "=", collapse = "; "),
    manifest        = sub(paste0("^", here::here(), "/?"), "", copy))
  f <- file.path(BASE_OUT_DIR, "run_log.csv")
  utils::write.table(row, f, sep = ",", row.names = FALSE, qmethod = "double",
                     append = file.exists(f), col.names = !file.exists(f))
  invisible(row)
}

stage_banner <- function(txt) {
  message("\n", strrep("=", 70), "\n  ", txt, "\n", strrep("=", 70))
}

# Column-name suffix for an extraction buffer radius, e.g. 0.02 -> "_2cm"
buffer_suffix <- function(b) paste0("_", round(b * 100), "cm")

# Exit status used by 03_models.R when it has just written the reference pool
# and stops for you to check it. run_all.R treats it as a pause, not an error.
STATUS_PAUSE <- 10L

# [NEW] End a stage early. From Rscript / run_all.R this exits with `status`;
# inside RStudio (interactive) it stops with a message instead, so it never
# tries to close RStudio.
end_stage <- function(msg, status = 0L) {
  message(msg)
  if (interactive()) {
    stop("Stage ended early (not an error): ", msg, call. = FALSE)
  }
  quit(save = "no", status = status)
}
