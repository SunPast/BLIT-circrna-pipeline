library(blit)

# -----------------------------
# Pipeline paths & config
# -----------------------------
pipeline_dir <- "/data/home/dingjia/BLIT-pipeline"

fqfile   <- file.path(pipeline_dir, "run_batch_from_qc/PHS003316/PHS003316.txt")
indir    <- "/data/home/dingjia/blit_fastp"
oudir    <- "/data/home/dingjia/blit_test/RNA_phs003316/result"
config   <- file.path(pipeline_dir, "config.sh")

# -----------------------------
# Threads per tool
# -----------------------------
nthreads_prog <- 16

# -----------------------------
# Step 0. Prepare
# -----------------------------
cat("[Step 0] Prepare reference & tools...\n")

# Indexes for CIRIquant
sink(file.path(oudir, "prepare_index.log"), split = TRUE)
exec("bash", file.path(pipeline_dir, "prepare/prepare_index.sh")) |>
  cmd_condaenv("CIRIquant") |>
  cmd_run()
sink()

# fastp: only if not exist
if (!dir.exists(indir)) {
  sink(file.path(oudir, "prepare_fastp.log"), split = TRUE)
  exec("bash", file.path(pipeline_dir, "prepare/prepare_fastp.sh")) |>
    cmd_condaenv("fastp") |>
    cmd_run()
  sink()
} else {
  cat("fastp output directory exists, skipping fastp prepare.\n")
}

# -----------------------------
# Step 1. Generate sample list
# -----------------------------
cat("[Step 1] Generate sample list...\n")
sink(file.path(oudir, "generate_sample_list.log"), split = TRUE)
exec(
  "python",
  file.path(pipeline_dir, "common/ll_fq.py"),
  indir,
  "--output", fqfile
) |>
  cmd_run()
sink()

samples <- readLines(fqfile)
cat("Total samples:", length(samples), "\n\n")

# -----------------------------
# Step 2. circRNA tools
# -----------------------------
tools <- list(
  CIRIquant = list(script = "CIRIquant/CIRIquant.sh", env = "CIRIquant"),
  FindCirc = list(script = "FindCirc/FindCirc.sh", env = "FindCirc"),
  Circexplorer2 = list(script = "Circexplorer2/Circexplorer2.sh", env = "Circexplorer2"),
  circRNA_finder = list(script = "circRNA_finder/circRNA_finder.sh", env = "circRNA_finder")
)

cat("[Step 2] Starting circRNA analysis (sequential mode)...\n")
cat("Tools:", paste(names(tools), collapse = ", "), "\n")
cat("Threads per tool:", nthreads_prog, "\n\n")

for (sample in samples) {
  cat("===> Processing sample:", sample, "\n")

  for (tool_name in names(tools)) {
    tool <- tools[[tool_name]]
    cat("----> Running", tool_name, "for", sample, "\n")

    log_file <- file.path(oudir, paste0(sample, ".", tool_name, ".log"))
    sink(log_file, split = TRUE)

    exec(
      "bash",
      file.path(pipeline_dir, tool$script),
      sample,
      indir,
      oudir,
      as.character(nthreads_prog),
      config
    ) |>
      cmd_condaenv(tool$env) |>
      cmd_run()

    sink()
  }

  cat("===> Sample finished:", sample, "\n\n")
}

cat("All circRNA tools completed (sequential mode)!\n")
cat("Now you can run the aggregation step manually.\n")
cat("Run: Rscript run_aggr.R\n")
