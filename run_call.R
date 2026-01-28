library(blit)
library(yaml)

generate_ciriquant_yaml <- function(cfg, outfile) {

  env_bin <- file.path(cfg$environments$env_CIRIquant, "bin")

  ciri_cfg <- list(
    name = "auto_build",
    tools = list(
      bwa       = file.path(env_bin, "bwa"),
      hisat2    = file.path(env_bin, "hisat2"),
      samtools  = file.path(env_bin, "samtools"),
      stringtie = file.path(env_bin, "stringtie")
    ),
    reference = list(
      fasta       = cfg$reference$fasta,
      gtf         = cfg$reference$gtf,
      bwa_index   = cfg$reference$fasta,
      hisat_index = cfg$reference$fasta
    )
  )

  yaml::write_yaml(ciri_cfg, outfile)
}

generate_config_sh <- function(cfg, outfile, ciri_yaml) {
  lines <- c(
    paste0("fasta=", cfg$reference$fasta),
    paste0("gtf=", cfg$reference$gtf),
    paste0("ann_ref=", cfg$reference$ann_ref),
    paste0("gdir=", cfg$reference$star_index),
    paste0("bt2_INDEX=", cfg$reference$bt2_INDEX),
    paste0("CIRI_config=", ciri_yaml)
  )
  writeLines(lines, outfile)
}

cfg <- yaml::read_yaml("input.yml")
ciri_yaml <- file.path(cfg$data$pipeline_dir, "CIRIquant/hg38.yml")
generate_ciriquant_yaml(cfg, ciri_yaml)
generate_config_sh(cfg, file.path(cfg$data$pipeline_dir, "config.sh"), ciri_yaml)

Sys.setenv(
  RAW_FASTQ_DIR = cfg$data$raw_fastq_dir,
  FASTP_OUTDIR = cfg$data$fastp_dir,
  REF_BASE = cfg$reference$base_dir,
  GENOME = cfg$reference$genome,
  REF_FASTA = cfg$reference$fasta,
  REF_GTF = cfg$reference$gtf,
  STAR_INDEX = cfg$reference$star_index,
  BT2_PREFIX = cfg$reference$bt2_INDEX,
  ENV_CIRCRNA_FINDER = cfg$environments$env_circRNA_finder,
  ENV_CIRIQUANT = cfg$environments$env_CIRIquant,
  ENV_FINDCIRC = cfg$environments$env_FindCirc,
  ENV_CE2 = cfg$environments$env_CE2,
  CIRI_CONFIG = ciri_yaml,
  THREADS = cfg$threads,
  MICROMAMBA = cfg$environments$micromamba
)

# -----------------------------
# Pipeline paths & config
# -----------------------------
pipeline_dir  <- cfg$data$pipeline_dir
fqfile        <- cfg$data$sample_list
indir         <- cfg$data$fastp_dir
outdir         <- cfg$data$result_dir
config        <- file.path(pipeline_dir, "config.sh")
nthreads_prog <- cfg$threads

# ----- log -----
run_call_log <- file.path(cfg$data$result_dir, "RUN_CALL.log")
dir.create(cfg$data$result_dir, showWarnings = FALSE, recursive = TRUE)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", args_full, value = TRUE)

if (Sys.getenv("BLIT_MASTER_LOG") == "" && length(script_arg) > 0) {
  log_file <- normalizePath(run_call_log, mustWork = FALSE)

  args_full <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args_full[grep("--file=", args_full)])
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

# -----------------------------
# Step 0. Prepare
# -----------------------------
cat("[Step 0] Prepare reference & tools...\n")

exec("bash", file.path(pipeline_dir, "prepare/prepare_index.sh")) |>
  cmd_run()

if (!dir.exists(indir)) {
  cat("fastp output directory not found, running fastp prepare...\n")
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

    exec(
      "bash",
      file.path(pipeline_dir, tool$script),
      sample,
      indir,
      outdir,
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
