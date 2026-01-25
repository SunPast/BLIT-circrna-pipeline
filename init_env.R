library(blit)

install_appmamba()

TUNA_CF <- "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge"
TUNA_BC <- "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda"
BASE <- "~/.local/share/R/blit/appmamba/envs"

# =======================
# QC: fastp
# =======================

appmamba(
  "create", "--yes",
  "-p", "~/.local/share/R/blit/appmamba/envs/fastp",
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "fastp"
)

# =======================
# FindCirc
# just:
# mamba create -n FindCirc -c conda-forge -c bioconda \
#   find_circ=1.2 bowtie2 samtools
# =======================

appmamba(
  "create", "--yes",
  "-p", file.path(BASE, "FindCirc"),
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "find_circ=1.2",
  "bowtie2",
  "samtools"
)

# =======================
# Circexplorer2
# just:
# mamba create -n Circexplorer2 -c conda-forge -c bioconda python=3.7 bwa
# pip install circexplorer2
# =======================

appmamba(
  "create", "--yes",
  "-p", file.path(BASE, "Circexplorer2"),
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "python=3.7",
  "bwa",
  "pip"
)

exec(
  file.path(BASE, "Circexplorer2/bin/pip"),
  "install",
  "circexplorer2"
) |>
  cmd_run()

# =======================
# CIRIquant
# just:
# mamba env create -f env.yml
# mamba install -n CIRIquant bwa=0.7.17 hisat2=2.2.0 \
#   stringtie=2.1.1 samtools=1.10
# mamba install -n CIRIquant r-base=3.6 r-optparse ...
# =======================

appmamba(
  "create", "--yes",
  "-p", file.path(BASE, "CIRIquant"),
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "python=2.7",
  "bwa=0.7.17",
  "hisat2=2.2.0",
  "stringtie=2.1.1",
  "samtools=1.10",
  "pip"
)

exec(
  file.path(BASE, "CIRIquant/bin/pip"),
  "install",
  "CIRIquant==1.1.2"
) |>
  cmd_run()

# =======================
# circRNA_finder
# just:
# mamba create -n circRNA_finder -c conda-forge -c bioconda \
#   circrna_finder star samtools
# =======================

appmamba(
  "create", "--yes",
  "-p", file.path(BASE, "circRNA_finder"),
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "circrna_finder",
  "star",
  "samtools"
)

cat("All environments installed (faithful to just).\n")

# =======================
# aggr environment (from just install)
# just:
# mamba install -c conda-forge -y r-base=4.2 r-data.table radian
# =======================

appmamba(
  "create", "--yes",
  "-p", file.path(BASE, "aggr"),
  "-c", TUNA_CF,
  "r-base=4.2",
  "r-data.table",
  "radian"
)

cat("aggr environment installed (R 4.2 + data.table + radian)\n")
