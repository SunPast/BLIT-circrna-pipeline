library(blit)

install_appmamba()

TUNA_CF <- "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge"
TUNA_BC <- "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda"

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
# Call: FindCirc
# (find_circ + bowtie2 + samtools)
# =======================

appmamba(
  "create", "--yes",
  "-p", "~/.local/share/R/blit/appmamba/envs/FindCirc",
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "bowtie2",
  "samtools"
)

# =======================
# Call: CIRIquant
# ⚠️ 必须 Python <= 3.8
# =======================

appmamba(
  "create", "--yes",
  "-p", "~/.local/share/R/blit/appmamba/envs/CIRIquant",
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "python=3.9",
  "ciriquant",
  "samtools",
  "bwa",
  "hisat2"
)

# =======================
# Call: CIRCexplorer2
# (pip only, requires Python)
# =======================

appmamba(
  "create", "--yes",
  "-p", "~/.local/share/R/blit/appmamba/envs/Circexplorer2",
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "python=3.9",
  "star",
  "bedtools",
  "samtools",
  "pip"
)

exec(
  "~/.local/share/R/blit/appmamba/envs/Circexplorer2/bin/pip",
  "install",
  "CIRCexplorer2"
) |>
  cmd_run()

# =======================
# Call: circRNA_finder
# =======================

appmamba(
  "create", "--yes",
  "-p", "~/.local/share/R/blit/appmamba/envs/circRNA_finder",
  "-c", TUNA_CF,
  "-c", TUNA_BC,
  "star",
  "perl",
  "samtools"
)

cat("All BLIT environments installed successfully.\n")
