#!/usr/bin/env bash
set -euo pipefail

REF_BASE=/data/home/dingjia/pipeline/
THREADS=20

LOG_FILE="${REF_BASE}/star_index_build.log"
MICROMAMBA="/data/home/dingjia/.local/bin/micromamba"

CIRCRNA_FINDER_ENV="/data/home/dingjia/.local/share/R/blit/appmamba/envs/circRNA_finder"
CIRIQUANT_ENV="/data/home/dingjia/.local/share/R/blit/appmamba/envs/CIRIquant"
FIND_CIRC_ENV="/data/home/dingjia/.local/share/R/blit/appmamba/envs/FindCirc"

FASTA=${REF_BASE}/GRCh38.primary_assembly.genome.fa
GTF=${REF_BASE}/gencode.v34.annotation.gtf

STAR_INDEX=${REF_BASE}/STAR_index_2.7.10b
BT2_PREFIX=${REF_BASE}/GRCh38.primary_assembly

echo "===> Genome index preparation"
echo "Log file: $LOG_FILE"

# bowtie2 (find_circ)
if [ ! -f "${BT2_PREFIX}.1.bt2" ]; then
  echo "===> bowtie2-build"
  ${MICROMAMBA} run -p "${FIND_CIRC_ENV}" \
    bowtie2-build \
    "${FASTA}" \
    "${BT2_PREFIX}"
else
  echo "===> Bowtie2 index exists at ${BT2_PREFIX}, skipping..."
fi

# bwa (CIRIquant)
if [ ! -f "${FASTA}.bwt" ]; then
  echo "===> bwa index"
  ${MICROMAMBA} run -p "${CIRIQUANT_ENV}" \
    bwa index -a bwtsw ${FASTA}
fi

# hisat2 (CIRIquant)
if [ ! -f "${FASTA}.1.ht2" ]; then
  echo "===> hisat2-build"
  ${MICROMAMBA} run -p "${CIRIQUANT_ENV}" \
    hisat2-build -p ${THREADS} ${FASTA} ${FASTA}
fi

# STAR (circRNA_finder)
if [ ! -d "${STAR_INDEX}" ] || [ ! -f "${STAR_INDEX}/Genome" ]; then
  echo "===> STAR genomeGenerate"
  echo "Start time: $(date)" | tee -a "$LOG_FILE"
  echo "FASTA: $FASTA" | tee -a "$LOG_FILE"
  echo "GTF: $GTF" | tee -a "$LOG_FILE"
  echo "Output: $STAR_INDEX" | tee -a "$LOG_FILE"
  echo "Threads: $THREADS" | tee -a "$LOG_FILE"
  echo "Environment: $CIRCRNA_FINDER_ENV" | tee -a "$LOG_FILE"

  rm -rf "${STAR_INDEX}"
  mkdir -p "${STAR_INDEX}"

  if [ ! -d "$CIRCRNA_FINDER_ENV" ]; then
    echo "Error: Environment path does not exist: $CIRCRNA_FINDER_ENV" | tee -a "$LOG_FILE"
    echo "Available environments:" | tee -a "$LOG_FILE"
    ${MICROMAMBA} env list | tee -a "$LOG_FILE"
    exit 1
  fi

  echo "Running STAR index construction..." | tee -a "$LOG_FILE"
  ${MICROMAMBA} run -p "$CIRCRNA_FINDER_ENV" \
    STAR \
      --runThreadN ${THREADS} \
      --runMode genomeGenerate \
      --genomeDir ${STAR_INDEX} \
      --genomeFastaFiles ${FASTA} \
      --sjdbGTFfile ${GTF} \
      --sjdbOverhang 100 \
      2>&1 | tee -a "$LOG_FILE"

  echo "End time: $(date)" | tee -a "$LOG_FILE"

  if [ -f "${STAR_INDEX}/Genome" ]; then
    echo "✓ STAR index built successfully" | tee -a "$LOG_FILE"
    echo "Index size: $(du -sh ${STAR_INDEX})" | tee -a "$LOG_FILE"
  else
    echo "✗ STAR index failed to build" | tee -a "$LOG_FILE"
    echo "Check log file: $LOG_FILE" | tee -a "$LOG_FILE"
    exit 1
  fi
else
  echo "===> STAR index exists at ${STAR_INDEX}, skipping..."
  echo "Index size: $(du -sh ${STAR_INDEX})"
fi

echo "Index preparation finished."
