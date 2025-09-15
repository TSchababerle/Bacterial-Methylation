#!/bin/bash
set -euo pipefail

# Usage check
if [[ $# -lt 3 ]]; then
  echo "Usage: bash methylation_pipeline.sh <assembly_dir> <methylation_bam_dir> <output_dir>"
  echo "  <assembly_dir> : directory with assemblies (.fasta)"
  echo "  <bam_dir>      : directory with methylation BAM files"
  echo "  <output_dir>   : directory for results"
  exit 1
fi

# Positional arguments
ASSEMBLY_DIR="$1"
BAM_DIR="$2"
OUTPUT_DIR="$3"

# Sub-directories
MAPPED_DIR="${OUTPUT_DIR}/Mapped_BAM"

# Create directories
mkdir -p "$OUTPUT_DIR" "$MAPPED_DIR"

echo "-----------------------------------"
echo "Methylation pipeline parameters:"
echo "  Assemblies : $ASSEMBLY_DIR"
echo "  BAM files  : $BAM_DIR"
echo "  Output dir : $OUTPUT_DIR"
echo "-----------------------------------"

# ID list
IDS=(
  "UT9728_MLa"
  "UT9728_MLb"
  "UT9728_IPa"
  "UT9728_IPb"
  "UT9728_TSa"
  "UT9728_TSb"
  "UT10237_JJa"
  "UT10237_JJb"
  "UT9728_OHa"
  "UT9728_OHb"
  "UT9728_HWa"
  "UT9728_HWb"
  "UT9728_MLc"
  "UT9728_IPc"
  "UT10237_TSc")

# Loop over each ID
for ID in "${IDS[@]}"; do
  ASSEMBLY="${ASSEMBLY_DIR}/${ID}.fasta"
  BAM_FILE="${BAM_DIR}/${ID}.bam"
  OUT_DIR="${OUTPUT_DIR}/${ID}"

  if [[ -f "$ASSEMBLY" && -f "$BAM_FILE" ]]; then
    echo "Processing $ID"
    echo "  Assembly: $ASSEMBLY"
    echo "  BAM file: $BAM_FILE"

    mkdir -p "$OUT_DIR"
    
	#load the conda environment
    eval "$(/path/to/conda/bin/conda shell.bash hook)"
    conda activate /path/to/conda/env
	
    #Map methylation BAM files to assembly
    samtools fastq "$BAM_FILE" -T MM,ML | \
       minimap2 -t 14 --secondary=no -ax map-ont -y "$ASSEMBLY" - | \
       samtools view -b | samtools sort -@ 10 -o "${MAPPED_DIR}/${ID}.mapped.bam"
    samtools index "${MAPPED_DIR}/${ID}.mapped.bam"
    BAM_FILE="${MAPPED_DIR}/${ID}.mapped.bam"
    
	conda deactivate
	
    ### MicrobeMod methylation calling
    conda activate microbemod
    MicrobeMod call_methylation \
      -b "$BAM_FILE" \
      -r "$ASSEMBLY" \
      -o "$OUT_DIR" \
      -t 10 

    #MicrobeMod annotate restriction modification systems pipeline
    MicrobeMod annotate_rm -f "$ASSEMBLY" -o "$OUT_DIR" -t 10

    conda deactivate

  else
    echo "Skipping $ID: Missing required files."
  fi
done
