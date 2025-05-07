#!/bin/bash

# Usage check
if [[ -z "$1" || -z "$2" || -z "$3" ]]; then
  echo "Usage: bash hybrid_pipeline.sh <long_reads_dir> <short_reads_dir> <output_dir>"
  exit 1
fi

# Input arguments
LR_FASTQ_DIR="$1"
SR_FASTQ_DIR="$2"
OUTPUT_DIR="$3"

#create output subdirectories
FASTQC_DIR="${OUTPUT_DIR}/1.fastqc"
FASTP_DIR="${OUTPUT_DIR}/2.fastp"
NANOPLOT_DIR="${OUTPUT_DIR}/3.nanoplot"
UNI_DIR="${OUTPUT_DIR}/4.unicycler"
QUAST_DIR="${OUTPUT_DIR}/5.quast"
BUSCO_DIR="${OUTPUT_DIR}/6.busco"
BWA_DIR="${OUTPUT_DIR}/7.bwa-mem"
POLYPOLISH_DIR="${OUTPUT_DIR}/8.polypolish"
BEROKKA_DIR="${OUTPUT_DIR}/9.berokka"
DNAAPLER_DIR="${OUTPUT_DIR}/10.dnaapler"
PROKKA_DIR="${OUTPUT_DIR}/11.prokka"

#create directories
mkdir -p $OUTPUT_DIR $FASTQC_DIR $NANOPLOT_DIR $FASTP_DIR $UNI_DIR \
         $POLYPOLISH_DIR $QUAST_DIR $BUSCO_DIR $PROKKA_DIR $BWA_DIR $BEROKKA_DIR $DNAAPLER_DIR

#define sample IDs
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
  SR_R1_FILE="${SR_FASTQ_DIR}/${ID}_R1.fastq.gz"
  SR_R2_FILE="${SR_FASTQ_DIR}/${ID}_R2.fastq.gz"
  LR_FILE="${LR_FASTQ_DIR}/${ID}.fastq.gz"

  if [[ -f "$SR_R1_FILE" && -f "$SR_R2_FILE" && -f "$LR_FILE" ]]; then
    echo "Processing $ID: Short reads ($SR_R1_FILE, $SR_R2_FILE), Long reads ($LR_FILE)"
	
    #activate general tools environment
    conda activate genome_assembly_pipeline
	
    #Run FastQC on short reads
    fastqc -t 4 -o "$FASTQC_DIR" "$SR_R1_FILE" "$SR_R2_FILE" &>> "${FASTQC_DIR}/fastqc_debug.log"

    #Run fastp for read trimming and filtering
    fastp -i "$SR_R1_FILE" -I "$SR_R2_FILE" -o "${FASTP_DIR}/${ID}_R1_trimmed.fastq.gz" -O "${FASTP_DIR}/${ID}_R2_trimmed.fastq.gz" --trim_poly_g &>> "${FASTP_DIR}/fastp_debug.log"

    #Run FastQC on trimmed reads
    fastqc -o "$FASTQC_DIR/${ID}" "${FASTP_DIR}/${ID}_R1_trimmed.fastq.gz" "${FASTP_DIR}/${ID}_R2_trimmed.fastq.gz" $>> "${FASTP_DIR}/fastp_debug.log"

    #Run NanoPlot on long reads
    NanoPlot --fastq "$LR_FILE" --outdir "${NANOPLOT_DIR}/${ID}" --threads 16 &>> "${NANOPLOT_DIR}/nanoplot_debug.log"

    #Run Unicycler on short reads and long reads
    unicycler -1 "${FASTP_DIR}/${ID}_R1_trimmed.fastq.gz" -2 "${FASTP_DIR}/${ID}_R2_trimmed.fastq.gz" -l "$LR_FILE" -o "${UNI_DIR}/${ID}"\
   	 --threads 16 \
	 --mode 'normal' \
	 --min_fasta_length '100' \
	 --linear_seqs '0' \
	 --min_kmer_frac '0.2' \
	 --max_kmer_frac '0.95' \
	 --kmer_count '10' \
	 --depth_filter '0.25' \
	 --start_gene_id '90.0' \
	 --start_gene_cov '95.0' \
	 --min_component_size '1000' \
	 --min_dead_end_size '1000' \
	 --scores '3,-6,-5,-2' \
	 --keep 1

    #Run Quast on raw assemblies
    quast.py -o "${QUAST_DIR}/${ID}_raw" "${UNI_DIR}/${ID}/assembly.fasta" &>> "${QUAST_DIR}/quast_debug.log"

    #Run Polypolish on assemblies
    bwa index "${UNI_DIR}/${ID}/assembly.fasta" -p "${BWA_DIR}/${ID}_indexed_assembly"
    bwa mem -t 16 -a "${BWA_DIR}/${ID}_indexed_assembly" "${FASTP_DIR}/${ID}_R1_trimmed.fastq.gz" > "${POLYPOLISH_DIR}/${ID}_alignments_1.sam"
    bwa mem -t 16 -a "${BWA_DIR}/${ID}_indexed_assembly" "${FASTP_DIR}/${ID}_R2_trimmed.fastq.gz" > "${POLYPOLISH_DIR}/${ID}_alignments_2.sam"
    polypolish filter --in1 "${POLYPOLISH_DIR}/${ID}_alignments_1.sam" --in2 "${POLYPOLISH_DIR}/${ID}_alignments_2.sam" --out1 "${POLYPOLISH_DIR}/${ID}_filtered_1.sam" --out2 "${POLYPOLISH_DIR}/${ID}_filtered_2.sam"
    polypolish polish "${UNI_DIR}/${ID}/assembly.fasta" "${POLYPOLISH_DIR}/${ID}_filtered_1.sam" "${POLYPOLISH_DIR}/${ID}_filtered_2.sam" > "${POLYPOLISH_DIR}/${ID}_polished.fasta"

    #Trim overlaps with berokka
    berokka --outdir "${BEROKKA_DIR}/${ID}" "${POLYPOLISH_DIR}/${ID}_polished.fasta" --force 

    #Rotate assemblies with dnaapler
    dnaapler all -i "${BEROKKA_DIR}/${ID}/02.trimmed.fa" -o "${DNAAPLER_DIR}/${ID}" -p "${ID}" -f -t 16

    #Run Quast on polished assemblies
    quast.py -o "${QUAST_DIR}/${ID}_polished" "${POLYPOLISH_DIR}/${ID}_polished.fasta" &>> "${QUAST_DIR}/quast_debug.log"

    #Annotate assembly with Prokka
    prokka --outdir "${PROKKA_DIR}/${ID}" --prefix "${ID}" "${POLYPOLISH_DIR}/${ID}_polished.fasta" &>> "${PROKKA_DIR}/prokka_debug.log"
    
    #Deactivate general_assembly_pipeline environment
    conda deactivate
	
    #Activate BUSCO environment
    conda activate busco_env
	
    #Run BUSCO on raw assemblies
    busco -i "${UNI_DIR}/${ID}/assembly.fasta" -o "${ID}_raw" --out_path "${BUSCO_DIR}/${ID}_raw" &>> "${BUSCO_DIR}/busco_debug.log" -l "/rsrch8/home/hlth_prof/taschababerle/sp25_paper/bacteria_odb10" -m genome -f --offline

    #Run BUSCO on polished assemblies
    busco -i "${POLYPOLISH_DIR}/${ID}_polished.fasta" -o "${ID}_polished" --out_path "${BUSCO_DIR}/${ID}_polished" &>> "${BUSCO_DIR}/busco_debug.log" -l "/rsrch8/home/hlth_prof/taschababerle/sp25_paper/bacteria_odb10" -m genome -f --offline
	
    #Deactivate busco environment
    conda deactivate
   else
     echo "skipping $ID: Missing required files."
   fi
done
