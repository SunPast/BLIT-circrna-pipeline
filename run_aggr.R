library(blit)
library(yaml)

cfg <- yaml::read_yaml("input.yml")

pipeline_dir  <- cfg$data$pipeline_dir
outdir        <- cfg$data$result_dir
fqfile       <- file.path(pipeline_dir, "PHS003316.txt")

aggr_dir <- file.path(outdir, "aggr")
dir.create(aggr_dir, showWarnings = FALSE, recursive = TRUE)

sample_anno <- fqfile

Sys.setenv(
  GTF_PATH  = cfg$reference$gtf,
  COMMON_PY = file.path(pipeline_dir, "common/common.py")
)

# ----- log -----
run_aggr_log <- file.path(outdir, "RUN_AGGR.log")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

args_full  <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", args_full, value = TRUE)

if (Sys.getenv("BLIT_MASTER_LOG") == "" && length(script_arg) > 0) {

  log_file <- normalizePath(run_aggr_log, mustWork = FALSE)

  script_path <- sub("--file=", "", script_arg)
  script_path <- normalizePath(script_path, mustWork = FALSE)

  cmd <- paste(
    shQuote(file.path(Sys.getenv("R_HOME"), "bin", "Rscript")),
    shQuote(script_path),
    ">>", shQuote(log_file), "2>&1"
  )

  Sys.setenv(BLIT_MASTER_LOG = "1")
  exec("bash", "-c", cmd) |> cmd_run()
  quit(save = "no")
}

# ----- check -----
cat("Checking if all tools have completed...\n")

check_completion <- function() {
  methods <- c("CIRI", "find_circ", "circexplorer2", "circRNA_finder")
  samples <- readLines(fqfile)

  all_complete <- TRUE
  for (sample in samples) {
    for (method in methods) {
      bed_file <- file.path(outdir, paste0(sample, ".", method, ".bed"))
      if (!file.exists(bed_file) || file.info(bed_file)$size == 0) {
        cat("Missing or empty:", basename(bed_file), "\n")
        all_complete <- FALSE
      }
    }
  }
  all_complete
}

if (!check_completion()) {
  stop("Not all tools finished. Fix missing BEDs first.")
}

cat("All tools completed. Starting aggregation...\n")

# ================================
# Step 1 — aggr_beds
# ================================
exec(
  "Rscript",
  file.path(pipeline_dir, "aggr/aggr_beds.R"),
  outdir,
  aggr_dir
) |>
  cmd_condaenv("aggr") |>
  cmd_run()

# ================================
# Step 2 — aggr_dataset
# ================================
exec(
  "Rscript",
  file.path(pipeline_dir, "aggr/aggr_dataset.R"),
  aggr_dir,
  aggr_dir,
  sample_anno
) |>
  cmd_condaenv("aggr") |>
  cmd_run()

cat("\nAggr finished\n")
cat("Aggregated results in:", aggr_dir, "\n")
