library(blit)

install_appmamba()

# ---------- QC ----------
appmamba("create", "--yes",
         "-p", "~/.local/share/R/blit/appmamba/envs/fastp",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda",
         "fastp")

# ---------- Call ----------
appmamba("create", "--yes",
         "-p", "~/.local/share/R/blit/appmamba/envs/FindCirc",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda",
         "bowtie2", "samtools")

appmamba("create", "--yes",
         "-p", "~/.local/share/R/blit/appmamba/envs/CIRIquant",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda",
         "bwa", "hisat2", "samtools")

appmamba("create", "--yes",
         "-p", "~/.local/share/R/blit/appmamba/envs/Circexplorer2",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda",
         "star", "bedtools", "python=3.10")

appmamba("create", "--yes",
         "-p", "~/.local/share/R/blit/appmamba/envs/circRNA_finder",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
         "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda",
         "star", "perl")
