#!/usr/bin/env bash
set -euo pipefail

# 0. Config
sample=$1
indir=$2
oudir=$3
ncpu=$4
prefix=${sample}

config=$5
source ${config}

ENV_PATH=${ENV_FindCirc}
MICROMAMBA=${MICROMAMBA}

outdir2=${oudir}/${prefix}.find_circ
mkdir -p ${outdir2}
cd ${outdir2}

echo "Start find_circ for ${sample} at $(date)"

if [ -f ../${prefix}.find_circ.bed ] && [ -s ../${prefix}.find_circ.bed ]; then
    echo "Final result exists and is not empty, skipped this sample."
    exit 0
fi

rm -f *.bam *.bai *.fastq *.fa *.log *.bed *.txt 2>/dev/null

export OMP_NUM_THREADS=${ncpu}
export OPENBLAS_NUM_THREADS=${ncpu}
export MKL_NUM_THREADS=${ncpu}

echo "1. Aligning reads..."

${MICROMAMBA} run -p ${ENV_PATH} bash -c "
set -euo pipefail

bowtie2 \
  -x ${bt2_INDEX} \
  -q -1 ${indir}/${sample}_1.fastq.gz -2 ${indir}/${sample}_2.fastq.gz \
  --threads ${ncpu} \
  --very-sensitive --score-min=C,-15,0 --reorder --mm \
  2> ${prefix}.bowtie2.log | \
samtools sort -@ ${ncpu} -o ${prefix}.bam -

samtools index -@ ${ncpu} ${prefix}.bam
"

echo "2. Fetching unmapped reads"

${MICROMAMBA} run -p ${ENV_PATH} bash -c "
set -euo pipefail

samtools view -@ ${ncpu} -hbf 4 ${prefix}.bam > ${prefix}.unmapped.bam
samtools index -@ ${ncpu} ${prefix}.unmapped.bam
"

echo "3. Splitting into anchors"

${MICROMAMBA} run -p ${ENV_PATH} unmapped2anchors.py \
    ${prefix}.unmapped.bam > ${prefix}.unmapped.fastq

[ -s ${prefix}.unmapped.fastq ] || { echo '[ERROR] unmapped.fastq empty'; exit 1; }

echo "4. Aligning anchors and running find_circ"

${MICROMAMBA} run -p ${ENV_PATH} bash -c "
set -euo pipefail

bowtie2 -q -U ${prefix}.unmapped.fastq -x ${bt2_INDEX} \
  --threads ${ncpu} \
  --reorder --mm --very-sensitive --score-min=C,-15,0 \
  2> ${prefix}.bowtie2.2nd.log | \
find_circ.py -G ${fasta} -n ${prefix} \
  --stats ${prefix}.sites.log \
  --reads ${prefix}.spliced_reads.fa \
  > ${prefix}.splice_sites.bed
"

[ -s ${prefix}.splice_sites.bed ] || { echo '[ERROR] splice_sites.bed not generated'; exit 1; }

grep CIRCULAR ${prefix}.splice_sites.bed | \
    grep -v chrM | \
    grep UNAMBIGUOUS_BP | \
    grep ANCHOR_UNIQUE | \
    maxlength.py 100000 \
    > ${prefix}.txt

awk -v OFS="\t" '{print $1,$2,$3,$6,$5,$4}' \
    ${prefix}.txt > ../${prefix}.find_circ.bed

rm -f *.bam *.bai *.fastq *.fa *.log *.bed *.txt 2>/dev/null

echo "Done for ${sample}, final result should be ${prefix}.find_circ.bed"
echo "End find_circ for ${sample} at $(date)"
