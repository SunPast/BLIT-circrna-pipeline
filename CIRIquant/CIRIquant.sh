#!/usr/bin/env bash
set -euo pipefail

sample=$1
indir=$2
oudir=$3
ncpu=$4
prefix=${sample}

config=$5
source ${config}

ENV_PATH=${ENV_CIRIQUANT}
MICROMAMBA=${MICROMAMBA}

outdir2=${oudir}/${prefix}.CIRI
mkdir -p ${outdir2}
cd ${outdir2}

echo "Start CIRIquant for ${sample} at $(date)"

if [ -f ../${prefix}.CIRI.bed ] && [ -s ../${prefix}.CIRI.bed ]; then
    echo "Final result exists and is not empty, skipped this sample."
    exit 0
fi

rm -rf align circ ${prefix}.gtf ${prefix}.log 2>/dev/null

export OMP_NUM_THREADS=${ncpu}
export OPENBLAS_NUM_THREADS=${ncpu}
export MKL_NUM_THREADS=${ncpu}

${MICROMAMBA} run -p ${ENV_PATH} \
    CIRIquant -t ${ncpu} \
        -1 ${indir}/${sample}_1.fastq.gz \
        -2 ${indir}/${sample}_2.fastq.gz \
        --config ${CIRI_config} \
        --no-gene \
        -o ${outdir2} \
        -p ${sample} \
        -v

if [ ! -s ${prefix}.gtf ]; then
    echo "[ERROR] CIRIquant failed: ${prefix}.gtf not generated" >&2
    exit 1
fi

grep -v "#" ${prefix}.gtf | awk '{print $14}' | cut -d '.' -f1 > ${prefix}.counts
grep -v "#" ${prefix}.gtf | \
    awk -v OFS="\t" '{gsub(/[";]/, "", $20); gsub(/[";]/, "", $22); print $1,$4-1,$5,$7,$20,$22}' \
    > ${prefix}.tmp

paste ${prefix}.tmp ${prefix}.counts > ../${prefix}.CIRI.bed

rm -f ${prefix}.tmp ${prefix}.counts ${prefix}.gtf ${prefix}.log 2>/dev/null
rm -rf align circ 2>/dev/null

echo "Done for ${sample}, final result should be ${prefix}.CIRI.bed"
echo "End CIRIquant for ${sample} at $(date)"
