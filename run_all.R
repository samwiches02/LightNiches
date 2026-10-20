# ══════════════════════════════════════════════════════════════════════════════
#  run_all.R — run the whole analysis, each stage in a FRESH R session
#
#  From a terminal in the project folder:
#     Rscript run_all.R          # all stages
#     Rscript run_all.R 03       # start at stage 03 (01-02 outputs must exist)
#     Rscript run_all.R 03 03    # run only stage 03
#
#  From RStudio: open run_all.R and click Source.
#
#  Why separate sessions: each stage can only see what earlier stages wrote to
#  output/. Nothing leaks from your workspace (master_script.R's
#  `if (exists("validation"))` blocks silently used whatever was left over).
#
#  PROVENANCE: [NEW] — written by Claude.
# ══════════════════════════════════════════════════════════════════════════════

STAGES <- c("01" = "scripts/01_inputs.R",
            "02" = "scripts/02_raytrace.R",
            "03" = "scripts/03_models.R",
            "04" = "scripts/04_figures.R",
            "07" = "scripts/07_taylor.R")
STATUS_PAUSE <- 10L

args  <- commandArgs(trailingOnly = TRUE)
from  <- if (length(args) >= 1) args[1] else names(STAGES)[1]
to    <- if (length(args) >= 2) args[2] else names(STAGES)[length(STAGES)]
run   <- STAGES[match(from, names(STAGES)):match(to, names(STAGES))]

root    <- here::here()
rscript <- file.path(R.home("bin"), "Rscript")

for (id in names(run)) {
  script <- file.path(root, run[[id]])
  message("\n>>> Stage ", id, ": ", run[[id]], "  (", format(Sys.time(), "%H:%M:%S"), ")")
  status <- system2(rscript, shQuote(script))
  if (status == STATUS_PAUSE) {
    message("\nPaused after stage ", id, ": check output/reference_pool.csv, ",
            "then run:  Rscript run_all.R ", id)
    break
  }
  if (status != 0) stop("Stage ", id, " failed (exit status ", status, "). ",
                        "Fix it, then resume with:  Rscript run_all.R ", id, call. = FALSE)
}
message("\nFinished.")
