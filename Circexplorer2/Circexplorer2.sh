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

ENV_PATH="/data/home/dingjia/.local/share/R/blit/appmamba/envs/Circexplorer2"
MICROMAMBA="/data/home/dingjia/.local/bin/micromamba"

#========================================================
outdir2=${oudir}/${prefix}.circexplorer2
mkdir -p ${outdir2}
cd ${outdir2}

echo "Start circexplorer2 for ${sample} at $(date)"

if [ -f ../${prefix}.circexplorer2.bed ] && [ -s ../${prefix}.circexplorer2.bed ]; then
    echo "Final result exists and is not empty, skipped this sample."
    exit 0
fi

rm -f *.sam *.txt *.log *.bed 2>/dev/null

export OMP_NUM_THREADS=${ncpu}
export OPENBLAS_NUM_THREADS=${ncpu}
export MKL_NUM_THREADS=${ncpu}

unmap_bwa_sam=${prefix}_unmapped_bwa.sam

echo "1. Aligning reads with BWA"

${MICROMAMBA} run -p ${ENV_PATH} bash -c "
set -euo pipefail

bwa mem -t ${ncpu} -T 19 \
  ${fasta} \
  ${indir}/${sample}_1.fastq.gz \
  ${indir}/${sample}_2.fastq.gz \
  > ${unmap_bwa_sam} 2> ${prefix}_bwa.log
"

[ -s ${unmap_bwa_sam} ] || { echo '[ERROR] BWA SAM not generated'; exit 1; }

echo "2. Running CIRCexplorer2 parse"

${MICROMAMBA} run -p ${ENV_PATH} bash -c "
set -euo pipefail

CIRCexplorer2 parse \
  -t BWA \
  -b ${prefix}_circ2_result.txt \
  ${unmap_bwa_sam} \
  > ${prefix}_parse.log
"

[ -s ${prefix}_circ2_result.txt ] || { echo '[ERROR] parse output missing'; exit 1; }

echo "3. Running CIRCexplorer2 annotate"

${MICROMAMBA} run -p ${ENV_PATH} bash -c "
set -euo pipefail

CIRCexplorer2 annotate \
  -r ${ann_ref} \
  -g ${fasta} \
  -b ${prefix}_circ2_result.txt \
  -o ${prefix}_circ2_result_ann.txt
"

[ -s ${prefix}_circ2_result_ann.txt ] || { echo '[ERROR] annotate output missing'; exit 1; }

awk -v OFS="\t" '{print $1,$2,$3,$6,$13}' \
  ${prefix}_circ2_result_ann.txt > ../${prefix}.circexplorer2.bed

rm -f *.sam *.txt *.log 2>/dev/null

echo "Done for ${sample}, final result should be ${prefix}.circexplorer2.bed"
echo "End circexplorer2 for ${sample} at $(date)"
