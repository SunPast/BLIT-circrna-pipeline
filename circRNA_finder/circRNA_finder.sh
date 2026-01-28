#!/usr/bin/env bash
#set -euxo pipefail

sample=$1
indir=$2
oudir=$3
ncpu=$4
prefix=${sample}

config=$5
source ${config}

ENV_PATH=${ENV_CIRCRNA_FINDER}
MICROMAMBA=${MICROMAMBA}

# NOTE:
# ncpu is controlled by the R pipeline (default = 4),
# so this script will NOT occupy excessive CPU resources.

outdir2=${oudir}/${prefix}.circRNA_finder
mkdir -p ${outdir2}
cd ${outdir2}

echo "Start circRNA_finder for ${sample} at $(date)"

if [ -f ../${prefix}.circRNA_finder.bed ] && [ -s ../${prefix}.circRNA_finder.bed ]; then
    echo "Final result exists and is not empty, skipped this sample."
    exit 0
elif [ -f ../${prefix}.circRNA_finder.bed ]; then
    echo "Final result file exists but is empty, re-run it."
else
    echo "Final result file does not exist, run it."
fi

# Clean working directory
rm -rf ./*

echo "1. Aligning reads with STAR..."

${MICROMAMBA} run -p ${ENV_PATH} \
  STAR \
  --readFilesIn ${indir}/${sample}_1.fastq.gz ${indir}/${sample}_2.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN ${ncpu} \
  --genomeDir ${gdir} \
  --chimSegmentMin 20 \
  --chimScoreMin 1 \
  --alignIntronMax 100000 \
  --chimOutType Junctions SeparateSAMold \
  --outFilterMismatchNmax 4 \
  --alignTranscriptsPerReadNmax 100000 \
  --outFilterMultimapNmax 2 \
  --outFileNamePrefix ${sample}. \
  > ${sample}.STAR.stdout.log 2> ${sample}.STAR.stderr.log

echo "2. Running circRNA_finder post-processing..."

${MICROMAMBA} run -p ${ENV_PATH} \
  postProcessStarAlignment.pl \
  --starDir ./ \
  --outDir ./ \
  > ${sample}.postprocess.log 2>&1

echo "3. Generating BED output..."

awk -v OFS="\t" '{print $1,$2,$3,$6,$5}' \
  ${prefix}.filteredJunctions.bed \
  > ../${prefix}.circRNA_finder.bed

rm -f \
  *.sam \
  *.bam \
  *.bai \
  *.junction \
  *.tab \
  Log.* \
  *.log \
  *_STARtmp \
  *_STARgenome \
  2>/dev/null

echo "Done for ${sample}, final result: ../${prefix}.circRNA_finder.bed"
echo "End circRNA_finder for ${sample} at $(date)"
