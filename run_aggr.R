library(blit)

pipeline_dir <- "/data/home/dingjia/BLIT-pipeline"
oudir    <- "/data/home/dingjia/blit_test/RNA_phs003316/result"
fqfile   <- file.path(pipeline_dir, "run_batch_from_qc/PHS003316/PHS003316.txt")

aggr_dir <- file.path(oudir, "aggr")
dir.create(aggr_dir, showWarnings = FALSE, recursive = TRUE)

sample_anno <- fqfile

# Set environment variables for R scripts
Sys.setenv(
  GTF_PATH = "/data/home/dingjia/pipeline/gencode.v34.annotation.gtf",
  COMMON_PY = "/data/home/dingjia/BLIT-pipeline/common/common.py"
)

# Use absolute path for Rscript
rscript_path <- "/opt/R-4.4.3/bin/Rscript"

# Verify Rscript exists
if (!file.exists(rscript_path)) {
  stop("Rscript not found at: ", rscript_path)
}

cat("Checking if all tools have completed...\n")

check_completion <- function() {
  methods <- c("CIRI", "find_circ", "circexplorer2", "circRNA_finder")
  samples <- readLines(fqfile)

  all_complete <- TRUE
  for (sample in samples) {
    for (method in methods) {
      bed_file <- file.path(oudir, paste0(sample, ".", method, ".bed"))
      if (!file.exists(bed_file) || file.info(bed_file)$size == 0) {
        cat("Missing or empty:", basename(bed_file), "\n")
        all_complete <- FALSE
      }
    }
  }
  return(all_complete)
}

if (check_completion()) {
  cat("All tools completed. Starting aggregation...\n")

  # 1) aggregate beds from 4 callers
  exec(
    rscript_path,
    file.path(pipeline_dir, "aggr/aggr_beds.R"),
    oudir,
    aggr_dir
  ) |>
    cmd_run()

  # 2) build aggregated dataset
  exec(
    rscript_path,
    file.path(pipeline_dir, "aggr/aggr_dataset.R"),
    aggr_dir,
    aggr_dir,
    sample_anno
  ) |>
    cmd_run()

  cat("BLIT pipeline completed successfully!\n")
  cat("Results saved to:", oudir, "\n")
  cat("Aggregated results in:", aggr_dir, "\n")
} else {
  cat("Not all tools have completed yet.\n")
  cat("Please wait for all .bed files to be generated before running aggregation.\n")
}
