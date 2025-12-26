#!/usr/bin/env bash
set -euo pipefail

REF_BASE=/home/dingjia/ref

THREADS=40

CIRCEXPLORER_ENV=Circexplorer2

FASTA=${REF_BASE}/GRCh38.primary_assembly.genome.fa
GTF=${REF_BASE}/gencode.v34.annotation.gtf
ANN_REF=${REF_BASE}/hg38_ref_all.txt

STAR_INDEX=${REF_BASE}/STAR_index_2.7.10b
BT2_PREFIX=${REF_BASE}/GRCh38.primary_assembly
CIRI_YML=${REF_BASE}/hg38.yml

mkdir -p ${REF_BASE}

echo "===> Step 1. Download genome fasta & gtf"

cd ${REF_BASE}

# genome fasta
if [ ! -f "${FASTA}" ]; then
  wget -O GRCh38.primary_assembly.genome.fa.gz \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
  gunzip GRCh38.primary_assembly.genome.fa.gz
fi

# gtf
if [ ! -f "${GTF}" ]; then
  wget -O gencode.v34.annotation.gtf.gz \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz
  gunzip gencode.v34.annotation.gtf.gz
fi

echo "===> Step 2. Circexplorer2 annotation"

if [ ! -f "${ANN_REF}" ]; then
  source activate ${CIRCEXPLORER_ENV}
  fetch_ucsc.py hg38 > ${ANN_REF}
  conda deactivate
fi

echo "===> Step 3. STAR index (circRNA_finder)"

if [ ! -d "${STAR_INDEX}" ]; then
  mkdir -p ${STAR_INDEX}
  STAR \
    --runThreadN ${THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${STAR_INDEX} \
    --genomeFastaFiles ${FASTA}
fi

echo "===> Step 4. bwa & hisat2 index (CIRIquant)"

if [ ! -f "${FASTA}.bwt" ]; then
  bwa index -a bwtsw ${FASTA}
fi

if [ ! -f "${FASTA}.1.ht2" ]; then
  hisat2-build -p ${THREADS} ${FASTA} ${FASTA}
fi

echo "===> Step 5. bowtie2 index (find_circ)"

if [ ! -f "${BT2_PREFIX}.1.bt2" ]; then
  bowtie2-build --threads ${THREADS} ${FASTA} ${BT2_PREFIX}
fi

echo "===> Step 6. generate CIRIquant hg38.yml template"

if [ ! -f "${CIRI_YML}" ]; then
cat << EOF > ${CIRI_YML}
genome:
  fasta: ${FASTA}
  hisat2_index: ${FASTA}
  bwa_index: ${FASTA}
annotation:
  gtf: ${GTF}
tools:
  bwa: bwa
  hisat2: hisat2
  samtools: samtools
EOF
fi

echo "All reference files and indexes are ready!"



RAW_DATA_DIR=/home/dingjia/blit_test/RNA_phs003316/raw
PROCESSED_DIR=/data/data2/IO_RNA/PHS003316/fq

mkdir -p ${PROCESSED_DIR}

for R1_FILE in ${RAW_DATA_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename ${R1_FILE} _1.fastq.gz)
    R1=${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz
    R2=${RAW_DATA_DIR}/${SAMPLE}_2.fastq.gz
    
    echo "Processing: ${SAMPLE}"
    
    fastp -i ${R1} -I ${R2} \
          -o ${PROCESSED_DIR}/${SAMPLE}_1.fastq.gz \
          -O ${PROCESSED_DIR}/${SAMPLE}_2.fastq.gz \
          -j ${PROCESSED_DIR}/${SAMPLE}_fastp.json \
          -h ${PROCESSED_DIR}/${SAMPLE}_fastp.html \
          --detect_adapter_for_pe \
          --cut_front --cut_tail \
          --cut_window_size 4 \
          --cut_mean_quality 20 \
          --length_required 50 \
          -w ${THREADS}
done

echo "Fastq preprocessing completed!"