#!/usr/bin/env bash
set -euo pipefail

REF_BASE=/data/home/dingjia/pipeline/
THREADS=40

FASTA=${REF_BASE}/GRCh38.primary_assembly.genome.fa
GTF=${REF_BASE}/gencode.v34.annotation.gtf

STAR_INDEX=${REF_BASE}/STAR_index_2.7.10b
BT2_PREFIX=${REF_BASE}/GRCh38.primary_assembly

echo "===> Genome index preparation"

# bowtie2 (find_circ)
if [ ! -f "${BT2_PREFIX}.1.bt2" ]; then
  echo "===> bowtie2-build"
  ~/.local/bin/micromamba run -n FindCirc \
    bowtie2-build \
    "${FASTA}" \
    "${BT2_PREFIX}"
else
  echo "===> Bowtie2 index exists at ${BT2_PREFIX}, skipping..."
fi

# bwa (CIRIquant)
if [ ! -f "${FASTA}.bwt" ]; then
  echo "===> bwa index"
  bwa index -a bwtsw ${FASTA}
fi

# hisat2 (CIRIquant)
if [ ! -f "${FASTA}.1.ht2" ]; then
  echo "===> hisat2-build"
  hisat2-build -p ${THREADS} ${FASTA} ${FASTA}
fi

# STAR (circRNA_finder)
if [ ! -d "${STAR_INDEX}" ]; then
  echo "===> STAR genomeGenerate"
  mkdir -p ${STAR_INDEX}
  STAR \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDEX} \
    --genomeFastaFiles ${FASTA}
fi

echo "Index preparation finished."
