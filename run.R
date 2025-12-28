library(blit)

pipeline_dir <- "/data/home/dingjia/BLIT-pipeline"

fqfile   <- file.path(pipeline_dir, "run_batch_from_qc/PHS003316/PHS003316.txt")
indir    <- "/data/home/dingjia/blit_fastp"
oudir    <- "/data/home/dingjia/blit_test/RNA_phs003316/result"
nthreads <- 20
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

nthreads_prog <- 8
njob <- max(1, nthreads %/% nthreads_prog)

tasks <- list()
for (sample in samples) {
  for (tool in names(tools)) {
    tasks <- c(
      tasks,
      list(
        exec(
          "bash",
          file.path(pipeline_dir, tools[[tool]]$script),
          sample,
          indir,
          oudir,
          as.character(nthreads_prog),
          config
        ) |>
          cmd_condaenv(tools[[tool]]$env)
      )
    )
  }
}

cmd_parallel(!!!tasks, threads = njob)

# ---------- Step 3. aggr ----------

aggr_dir <- file.path(oudir, "aggr")
dir.create(aggr_dir, showWarnings = FALSE, recursive = TRUE)

sample_anno <- fqfile

# 1) aggregate beds from 4 callers
exec(
  "Rscript",
  file.path(pipeline_dir, "aggr/aggr_beds.R"),
  oudir,
  aggr_dir
) |>
  cmd_run()

# 2) build aggregated dataset
exec(
  "Rscript",
  file.path(pipeline_dir, "aggr/aggr_dataset.R"),
  aggr_dir,
  aggr_dir,
  sample_anno
) |>
  cmd_run()
