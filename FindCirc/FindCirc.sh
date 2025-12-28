#!/usr/bin/env bash
#set -euxo pipefail

# 0. Config
sample=$1
indir=$2
oudir=$3
ncpu=$4
prefix=${sample}

config=$5
source ${config}

ENV_PATH="/data/home/dingjia/.local/share/R/blit/appmamba/envs/FindCirc"
MICROMAMBA="/data/home/dingjia/.local/bin/micromamba"

#========================================================
outdir2=${oudir}/${prefix}.find_circ
mkdir -p ${outdir2}
cd ${outdir2}

echo "Start find_circ for ${sample} at `date`"

if [ -f ../${prefix}.find_circ.bed ] && [ -s ../${prefix}.find_circ.bed ]; then
    echo "Final result exists and is not empty, skipped this sample."
    exit 0
elif [ -f ../${prefix}.find_circ.bed ]; then
    echo "Final result file exists but is empty, re-run it."
else
    echo "Final result file does not exist, run it."
fi

rm -rf ${outdir2}/*

echo "1. Aligning reads..."

${MICROMAMBA} run -p ${ENV_PATH} \
    bowtie2 \
    -x ${bt2_INDEX} \
    -q -1 ${indir}/${sample}_1.fastq.gz -2 ${indir}/${sample}_2.fastq.gz \
    --threads ${ncpu} \
    --very-sensitive --score-min=C,-15,0 --reorder --mm \
    2> ${prefix}.bowtie2.log | \
    ${MICROMAMBA} run -p ${ENV_PATH} samtools sort -@ ${ncpu} -o ${prefix}.bam -

${MICROMAMBA} run -p ${ENV_PATH} samtools index -@ ${ncpu} ${prefix}.bam

echo "2. Fetching the unmapped reads"
${MICROMAMBA} run -p ${ENV_PATH} samtools view -@ ${ncpu} -hbf 4 ${prefix}.bam > ${prefix}.unmapped.bam
${MICROMAMBA} run -p ${ENV_PATH} samtools index -@ ${ncpu} ${prefix}.unmapped.bam

echo "3. Splitting into anchors"
${MICROMAMBA} run -p ${ENV_PATH} unmapped2anchors.py ${prefix}.unmapped.bam > ${prefix}.unmapped.fastq

echo "4. Aligning anchors and piping through find_circ"

${MICROMAMBA} run -p ${ENV_PATH} bowtie2 -q -U ${prefix}.unmapped.fastq -x ${bt2_INDEX} --threads ${ncpu} \
    --reorder --mm --very-sensitive --score-min=C,-15,0 2> ${prefix}.bowtie2.2nd.log | \
    ${MICROMAMBA} run -p ${ENV_PATH} find_circ.py -G ${fasta} -n ${prefix} \
    --stats ${prefix}.sites.log \
    --reads ${prefix}.spliced_reads.fa \
    > ${prefix}.splice_sites.bed

${MICROMAMBA} run -p ${ENV_PATH} grep CIRCULAR ${prefix}.splice_sites.bed | \
    ${MICROMAMBA} run -p ${ENV_PATH} grep -v chrM | \
    ${MICROMAMBA} run -p ${ENV_PATH} grep UNAMBIGUOUS_BP | ${MICROMAMBA} run -p ${ENV_PATH} grep ANCHOR_UNIQUE | \
    ${MICROMAMBA} run -p ${ENV_PATH} maxlength.py 100000 \
    > ${prefix}.txt

${MICROMAMBA} run -p ${ENV_PATH} awk -v OFS="\t" '{print $1,$2,$3,$6,$5,$4}' ${prefix}.txt > ../${prefix}.find_circ.bed

rm -f *.bam *.bai *.fastq *.fa *.log *.bed *.txt 2>/dev/null

echo "Done for ${sample}, final result should be ${prefix}.find_circ.bed"
echo "End find_circ for ${sample} at `date`"
