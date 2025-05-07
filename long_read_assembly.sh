#!/bin/bash

# Usage check
if [[ -z "$1" || -z "$2" ]]; then
  echo "Usage: bash hybrid_pipeline.sh <long_reads_dir> <output_dir>"
  exit 1
fi

# Input arguments
LR_FASTQ_DIR="$1"
OUTPUT_DIR="$2"

# Create output subdirectories
FASTQC_DIR="${OUTPUT_DIR}/1.fastqc"
NANOPLOT_DIR="${OUTPUT_DIR}/2.nanoplot"
FLYE_DIR="${OUTPUT_DIR}/3.flye"
MEDAKA_DIR="${OUTPUT_DIR}/6.medaka"
QUAST_DIR="${OUTPUT_DIR}/4.quast"
BUSCO_DIR="${OUTPUT_DIR}/5.busco"
PROKKA_DIR="${OUTPUT_DIR}/9.prokka"
BEROKKA_DIR="${OUTPUT_DIR}/7.berokka"
DNAAPLER_DIR="${OUTPUT_DIR}/8.dnaapler"

mkdir -p "$OUTPUT_DIR" "$FASTQC_DIR" "$NANOPLOT_DIR" "$FLYE_DIR" "$MEDAKA_DIR" "$QUAST_DIR" \
         "$BUSCO_DIR" "$PROKKA_DIR" "$BEROKKA_DIR" "$DNAAPLER_DIR"

# Define sample IDs
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
  "UT10237_TSc"
)

# Load conda base once
eval "$(/path/to/miniforge3/bin/conda shell.bash hook)"

# Process each identifier
for ID in "${IDS[@]}"; do
  LR_FILE="${LR_FASTQ_DIR}/${ID}.fastq.gz"

  if [[ -f "$LR_FILE" ]]; then
    echo "Processing $ID: Long reads ($LR_FILE)"

    # Activate general tools environment
    conda activate genome_assembly_pipeline

    # Run FastQC on long reads
    fastqc -t 4 -o "$FASTQC_DIR" "$LR_FILE" &>> "${FASTQC_DIR}/fastqc_debug.log"

    # Run NanoPlot on long reads
    NanoPlot --fastq "$LR_FILE" --outdir "${NANOPLOT_DIR}/${ID}" --threads 16 &>> "${NANOPLOT_DIR}/nanoplot_debug.log"

    # Assemble long reads with Flye
    flye --nano-hq "$LR_FILE" --scaffold -g 2.5m --out-dir "${FLYE_DIR}/${ID}" --threads 16 &>> "${FLYE_DIR}/flye_debug.log"

    # Run Quast on raw assemblies
    quast.py -o "${QUAST_DIR}/${ID}_raw" "${FLYE_DIR}/${ID}/assembly.fasta" &>> "${QUAST_DIR}/quast_debug.log"

    # Polish assembly using Medaka
    medaka_consensus -i "$LR_FILE" -d "${FLYE_DIR}/${ID}/assembly.fasta" -o "${MEDAKA_DIR}/${ID}" -t 16 -f &>> "${MEDAKA_DIR}/medaka_debug.log"

    # Trim overlaps with Berokka
    berokka --outdir "${BEROKKA_DIR}/${ID}" "${MEDAKA_DIR}/${ID}/consensus.fasta" --force &>> "${BEROKKA_DIR}/berokka_debug.log"

    # Rotate assembly using dnaapler
    dnaapler all -i "${BEROKKA_DIR}/${ID}/02.trimmed.fa" -o "${DNAAPLER_DIR}/${ID}" -p "${ID}" -f -t 16

    # Run Quast on polished assemblies
    quast.py -o "${QUAST_DIR}/${ID}_polished" "${DNAAPLER_DIR}/${ID}/${ID}_reoriented.fasta" &>> "${QUAST_DIR}/quast_debug.log"

    # Run Prokka to annotate assembly
    prokka --outdir "${PROKKA_DIR}/${ID}" --prefix "${ID}" "${DNAAPLER_DIR}/${ID}/${ID}_reoriented.fasta" &>> "${PROKKA_DIR}/prokka_debug.log"

    # Deactivate general tools environment
    conda deactivate

    # Activate BUSCO environment
    conda activate busco_env

    # Run BUSCO on raw assemblies
    busco -i "${FLYE_DIR}/${ID}/assembly.fasta" -o "${ID}_raw" --out_path "${BUSCO_DIR}/${ID}_raw" -l "/path/to/bacteria_odb10" -m genome -f --offline &>> "${BUSCO_DIR}/busco_debug.log"

    # Run BUSCO on polished assemblies
    busco -i "${DNAAPLER_DIR}/${ID}/${ID}_reoriented.fasta"  -o "${ID}_polished" --out_path "${BUSCO_DIR}/${ID}_polished" -l "/path/to/bacteria_odb10" -m genome -f --offline &>> "${BUSCO_DIR}/busco_debug.log"

    # Deactivate BUSCO environment
    conda deactivate
  else
    echo "Skipping $ID: Missing required files."
  fi
done
