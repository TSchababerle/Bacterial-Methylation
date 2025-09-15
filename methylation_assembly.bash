#!/bin/bash

# Usage check
if [[ $# -lt 4 ]]; then
  echo "Usage: bash hybrid_pipeline.sh <long_reads_dir> <short_reads_dir> <methylation_bam_dir> <output_dir>"
  exit 1
fi

# Positional arguments
LR_DIR="$1"        # Long reads
SR_DIR="$2"        # Short reads
BAM_DIR="$3"       # Pre-generated Methylation BAM files
OUTPUT_DIR="$4"    # Output directory

# Make sure output dir exists
mkdir -p "$OUTPUT_DIR"

echo "-----------------------------------"
echo "Hybrid pipeline parameters:"
echo "  Long reads : $LR_DIR"
echo "  Short reads: $SR_DIR"
echo "  Methylation BAM files: $BAM_DIR"
echo "  Output dir : $OUTPUT_DIR"
echo "-----------------------------------"

# Sub-directories
ASSEMBLY_DIR="${OUTPUT_DIR}/Raw_Assemblies"
POLISHED_DIR="${OUTPUT_DIR}/Polished_Assemblies"
MAPPED_DIR="${OUTPUT_DIR}/Mapped_BAM"
BWA_DIR="${OUTPUT_DIR}/BWA_Indexed_Genomes"
MICROBEMOD_DIR="${OUTPUT_DIR}/Microbemod"
PROKKA_DIR="${OUTPUT_DIR}/Prokka"

# Create directories
mkdir -p "$ASSEMBLY_DIR" "$POLISHED_DIR" "$MAPPED_DIR" "$BWA_DIR" "$MICROBEMOD_DIR" "$PROKKA_DIR"

# Define sample IDs (edit this list for your dataset)
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

# Loop over each sample
for ID in "${IDS[@]}"; do
  SR_R1_FILE="${SR_DIR}/${ID}_R1.fastq.gz"
  SR_R2_FILE="${SR_DIR}/${ID}_R2.fastq.gz"
  LR_FILE="${LR_DIR}/${ID}.fastq.gz"

  if [[ -f "$SR_R1_FILE" && -f "$SR_R2_FILE" && -f "$LR_FILE" ]]; then
    echo "Processing $ID"

    ### Flye assembly
    conda activate flye
    flye --nano-hq "$LR_FILE" --scaffold -g 2.5m --out-dir "${ASSEMBLY_DIR}/${ID}" --threads 16
    conda deactivate

    ### Polypolish (requires bwa + samtools)
    conda activate genome_assembly_pipeline
    bwa index "${ASSEMBLY_DIR}/${ID}/assembly.fasta" -p "${BWA_DIR}/${ID}_indexed_assembly"
    bwa mem -t 16 -a "${BWA_DIR}/${ID}_indexed_assembly" "$SR_R1_FILE" > "${POLISHED_DIR}/${ID}_alignments_1.sam"
    bwa mem -t 16 -a "${BWA_DIR}/${ID}_indexed_assembly" "$SR_R2_FILE" > "${POLISHED_DIR}/${ID}_alignments_2.sam"
    polypolish filter \
      --in1 "${POLISHED_DIR}/${ID}_alignments_1.sam" \
      --in2 "${POLISHED_DIR}/${ID}_alignments_2.sam" \
      --out1 "${POLISHED_DIR}/${ID}_filtered_1.sam" \
      --out2 "${POLISHED_DIR}/${ID}_filtered_2.sam"
    polypolish polish \
      "${ASSEMBLY_DIR}/${ID}/assembly.fasta" \
      "${POLISHED_DIR}/${ID}_filtered_1.sam" \
      "${POLISHED_DIR}/${ID}_filtered_2.sam" \
      > "${POLISHED_DIR}/${ID}_polished.fasta"

    ### Minimap2 mapping
    samtools fastq "${BAM_DIR}/${ID}.bam" -T MM,ML | \
      minimap2 -t 14 --secondary=no -ax map-ont -y "${POLISHED_DIR}/${ID}_polished.fasta" - | \
      samtools view -b | samtools sort -@ 10 -o "${MAPPED_DIR}/${ID}.mapped.bam"
    samtools index "${MAPPED_DIR}/${ID}.mapped.bam"
    conda deactivate

    ### Prokka annotation
    conda activate prokka
    prokka --outdir "${PROKKA_DIR}/${ID}" --prefix "${ID}" "${POLISHED_DIR}/${ID}_polished.fasta"
    conda deactivate

    ### MicrobeMod methylation calling
    conda activate microbemod
    MicrobeMod call_methylation \
      -b "${MAPPED_DIR}/${ID}.mapped.bam" \
      -r "${POLISHED_DIR}/${ID}_polished.fasta" \
      -o "${MICROBEMOD_DIR}/${ID}" -t 10

    MicrobeMod annotate_rm \
      -f "${POLISHED_DIR}/${ID}_polished.fasta" \
      -o "${MICROBEMOD_DIR}/${ID}" -t 10
    conda deactivate

  else
    echo "Skipping $ID: Missing required files."
  fi
done
