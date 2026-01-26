library(blit)

# -----------------------------
# Pipeline paths & config
# -----------------------------
pipeline_dir <- "/data/home/dingjia/BLIT-pipeline"

fqfile   <- file.path(pipeline_dir, "PHS003316.txt")
indir    <- "/data/home/dingjia/blit_fastp"
oudir    <- "/data/home/dingjia/blit_test/RNA_phs003316/result"
config   <- file.path(pipeline_dir, "config.sh")

nthreads_prog <- 48

# -----------------------------
# Step 0. Prepare
# -----------------------------
cat("[Step 0] Prepare reference & tools...\n")

exec("bash", file.path(pipeline_dir, "prepare/prepare_index.sh")) |>
  cmd_condaenv("CIRIquant") |>
  cmd_run()

if (!dir.exists(indir)) {
  exec("bash", file.path(pipeline_dir, "prepare/prepare_fastp.sh")) |>
    cmd_condaenv("fastp") |>
    cmd_run()
} else {
  cat("fastp output directory exists, skipping fastp prepare.\n")
}

# -----------------------------
# Step 1. Generate sample list
# -----------------------------
cat("[Step 1] Generate sample list...\n")
exec(
  "python",
  file.path(pipeline_dir, "common/ll_fq.py"),
  indir,
  "--output", fqfile
) |>
  cmd_run()

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

  }

  cat("===> Sample finished:", sample, "\n\n")
}

cat("All circRNA tools completed (sequential mode)!\n")
cat("Now you can run the aggregation step manually.\n")
cat("please run: run_aggr.R\n")
