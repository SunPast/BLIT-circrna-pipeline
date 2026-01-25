# CircRNA Identification and Integration Pipeline (BLIT Version)

## Purpose

This project is a **BLIT-integrated refactor** of the [original circRNA pipeline](https://github.com/OncoHarmony-Network/circrna-pipeline). The biological workflow is unchanged, but:

- All conda/mamba environment creation is now handled by BLIT (`appmamba`)
- All shell orchestration is replaced with BLIT functions like `exec()` and `cmd_condaenv()`
- The pipeline demonstrates how BLIT can manage a real multi-tool bioinformatics workflow

The pipeline detects circRNAs from paired-end FASTQ files using four methods:

- CIRIquant
- Circexplorer2
- find_circ
- circRNA_finder

Results from all methods are aggregated using an ensemble approach.

------

## ⚠️ Path Localization Required

This repository contains hard-coded paths from the author's machine. Before running, modify the following paths according to your local environment:

| File                               | PATHs to Modify                                            |
| :--------------------------------- | :--------------------------------------------------------- |
| `prepare/prepare_fastp.sh`         | Fastq input/output paths                                   |
| `prepare/prepare_index.sh`         | Genome, index, and environment paths                       |
| `Circexplorer2/Circexplorer2.sh`   | Environment paths                                          |
| `circRNA_finder/circRNA_finder.sh` | Environment paths                                          |
| `CIRIquant/CIRIquant.sh`           | Environment paths                                          |
| `CIRIquant/hg38.yml`               | Reference paths and tool PATHs in environments             |
| `FindCirc/FindCirc.sh`             | Environment paths                                          |
| `run_call.R`                       | sample list, fastp-cleaned files, output path, config file |
| `run_aggr.R`                       | pipeline directory, sample list, output path               |

All appmamba environments are created at: `~/.local/share/R/blit/appmamba/envs/`

------

## Step 1. Install Required Environments

Run in R:

```
source("init_env.R")
```

This reproduces the original justfile installations using BLIT.

Environments created:

| Environment    | Software                                                     |
| :------------- | :----------------------------------------------------------- |
| fastp          | fastp                                                        |
| FindCirc       | find_circ=1.2, bowtie2, samtools                             |
| Circexplorer2  | python=3.7, bwa, pip, circexplorer2                          |
| CIRIquant      | python=2.7, bwa=0.7.17, hisat2=2.2.0, stringtie=2.1.1, samtools=1.10, CIRIquant=1.1.2 |
| circRNA_finder | circrna_finder, star, samtools                               |
| aggr           | r-base=4.2, r-data.table, radian                             |

------

## Step 2. Prepare References and Configuration

1. **Prepare genome files:**

   - Genome FASTA (e.g., `GRCh38.primary_assembly.genome.fa`)
   - Annotation GTF (e.g., `gencode.v34.annotation.gtf`)

2. **For Circexplorer2:**
   Download `hg38_ref_all.txt` using `fetch_ucsc.py` inside its environment.

3. **Generate alignment indexes:**

   ```bash
   bash prepare/prepare_index.sh
   ```

4. **Configure paths in:**

   - `config.sh` - Main configuration
   - `CIRIquant/hg38.yml` - CIRIquant-specific settings

------

## Step 3. Preprocess FASTQ Files

Run quality control and adapter trimming with fastp:

```bash
bash prepare/prepare_fastp.sh
```

FASTQ naming convention:

```
*_1.fastq.gz  # Read 1
*_2.fastq.gz  # Read 2
```

------

## Step 4. Run circRNA Detection

Run in R:

```R
source("run_call.R")
```

This script:

1. Generates sample list from processed FASTQ files
2. Runs all four circRNA callers sequentially
3. Each tool runs in its BLIT-managed environment

Output per sample:

```
sample.CIRI.bed
sample.circexplorer2.bed
sample.find_circ.bed
sample.circRNA_finder.bed
```

------

## Step 5. Aggregate Results

Run in R:

```R
source("run_aggr.R")
```

This step:

1. Aggregates results from all four callers
2. Annotates circRNAs using GTF
3. Produces final integrated dataset

------

## Important Notes

1. **CIRIquant Installation:** This legacy tool only supports Python 2. It is installed via `pip install CIRIquant==1.1.2` inside a Python 2 environment.
2. **Environment Location:** All conda environments are installed at `~/.local/share/R/blit/appmamba/envs/`
3. **Resource Requirements:** Ensure sufficient disk space (approximately 20GB for environments and indexes)

------

## References

- Original pipeline: [OncoHarmony-Network/circrna-pipeline](https://github.com/OncoHarmony-Network/circrna-pipeline)
- Pipeline inspiration: [Yelab2020/ICBcircSig](https://github.com/Yelab2020/ICBcircSig), [nf-core/circrna](https://github.com/nf-core/circrna)
- BLIT: [WangLabCSU/BLIT](https://github.com/WangLabCSU/BLIT)
