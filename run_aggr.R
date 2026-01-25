library(blit)

pipeline_dir <- "/data/home/dingjia/BLIT-pipeline"
outdir        <- "/data/home/dingjia/blit_test/RNA_phs003316/result"
fqfile       <- file.path(pipeline_dir, "PHS003316.txt")

aggr_dir <- file.path(outdir, "aggr")
dir.create(aggr_dir, showWarnings = FALSE, recursive = TRUE)

sample_anno <- fqfile

Sys.setenv(
  GTF_PATH  = "/data/home/dingjia/pipeline/gencode.v34.annotation.gtf",
  COMMON_PY = "/data/home/dingjia/BLIT-pipeline/common/common.py"
)

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
