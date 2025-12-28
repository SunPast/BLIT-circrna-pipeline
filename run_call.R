library(blit)

pipeline_dir <- "/data/home/dingjia/BLIT-pipeline"

fqfile   <- file.path(pipeline_dir, "run_batch_from_qc/PHS003316/PHS003316.txt")
indir    <- "/data/home/dingjia/blit_fastp"
oudir    <- "/data/home/dingjia/blit_test/RNA_phs003316/result"
nthreads <- 40
config   <- file.path(pipeline_dir, "config.sh")

# ---------- Step 0. reference & QC ----------
exec("bash", file.path(pipeline_dir, "prepare/prepare_circexplorer.sh")) |>
  cmd_condaenv("Circexplorer2") |>
  cmd_run()

exec("bash", file.path(pipeline_dir, "prepare/prepare_index.sh")) |>
  cmd_condaenv("CIRIquant") |>
  cmd_run()

exec("bash", file.path(pipeline_dir, "prepare/prepare_fastp.sh")) |>
  cmd_condaenv("fastp") |>
  cmd_run()

# ---------- Step 1. sample list ----------
exec(
  "python",
  file.path(pipeline_dir, "common/ll_fq.py"),
  indir,
  "--output", fqfile
) |>
  cmd_run()

samples <- readLines(fqfile)

# ---------- Step 2. circRNA tools ----------
tools <- list(
  CIRIquant = list(script = "CIRIquant/CIRIquant.sh", env = "CIRIquant"),
  FindCirc = list(script = "FindCirc/FindCirc.sh", env = "FindCirc"),
  Circexplorer2 = list(script = "Circexplorer2/Circexplorer2.sh", env = "Circexplorer2"),
  circRNA_finder = list(script = "circRNA_finder/circRNA_finder.sh", env = "circRNA_finder")
)

nthreads_prog <- 4

cat("Starting circRNA analysis (sequential mode)...\n")
cat("Samples:", length(samples), "\n")
cat("Tools:", paste(names(tools), collapse = ", "), "\n")
cat("Threads per tool:", nthreads_prog, "\n\n")

for (sample in samples) {
  cat("===> Processing sample:", sample, "\n")

  for (tool in names(tools)) {
    cat("----> Running", tool, "for", sample, "\n")

    exec(
      "bash",
      file.path(pipeline_dir, tools[[tool]]$script),
      sample,
      indir,
      oudir,
      as.character(nthreads_prog),
      config
    ) |>
      cmd_condaenv(tools[[tool]]$env) |>
      cmd_run()
  }

  cat("===> Sample finished:", sample, "\n\n")
}

cat("All circRNA tools completed (sequential mode)!\n")
cat("Now you can run the aggregation step manually.\n")
cat("Run: Rscript run_aggr.R\n")
