#!/bin/bash
set -euo pipefail

RAW_DIR=${RAW_FASTQ_DIR}
OUT_DIR=${FASTP_DIR}

mkdir -p "$OUT_DIR"

cd "$RAW_DIR" || { echo "Failed to change directory to $RAW_DIR"; exit 1; }

echo "Starting fastp processing..."

# Find R1 files and process
find . -type f \( -name "*_1.fastq.gz" -o -name "*_R1.fastq.gz" \) | while read -r r1; do
    # Remove leading ./
    r1="${r1#./}"

    # Determine naming pattern and set r2 filename
    if [[ "$r1" == *_1.fastq.gz ]]; then
        r2="${r1/_1.fastq.gz/_2.fastq.gz}"
        prefix=$(basename "$r1" _1.fastq.gz)
        output_r1="${prefix}_1.fastq.gz"
        output_r2="${prefix}_2.fastq.gz"
    elif [[ "$r1" == *_R1.fastq.gz ]]; then
        r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
        prefix=$(basename "$r1" _R1.fastq.gz)
        output_r1="${prefix}_R1.fastq.gz"
        output_r2="${prefix}_R2.fastq.gz"
    else
        echo "Skipping unexpected file: $r1"
        continue
    fi

    # Check if output already exists
    if [[ -f "$OUT_DIR/$output_r1" ]] && [[ -f "$OUT_DIR/$output_r2" ]]; then
        echo "Output already exists for $prefix, skipping..."
        continue
    fi

    # Check if R2 file exists
    if [ -f "$r2" ]; then
        echo "Processing paired-end: $r1 and $r2"

        base_r1=$(basename "$r1")
        base_r2=$(basename "$r2")

        # Run fastp
        fastp -i "$r1" -I "$r2" \
              -o "$OUT_DIR/$base_r1" \
              -O "$OUT_DIR/$base_r2" \
              -w 40 \
              -j "$OUT_DIR/${prefix}.json" \
              -h "$OUT_DIR/${prefix}.html" \
              --detect_adapter_for_pe

        # Verify output
        if [[ -f "$OUT_DIR/$base_r1" ]] && [[ -f "$OUT_DIR/$base_r2" ]]; then
            echo "Successfully processed: $prefix"
        else
            echo "Warning: Output files for $prefix may not have been created properly"
        fi
    else
        echo "Error: Paired file $r2 not found for $r1"
        # Optionally process as single-end or skip
        # fastp -i "$r1" -o "$OUT_DIR/$base_r1" ...
    fi
done

echo "Processing completed. Output directory: $OUT_DIR"
