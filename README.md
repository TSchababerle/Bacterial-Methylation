# Bacterial-Methylation
This repository contains two bash scripts for bacterial genome assembly and a script for parsing genome features and annotating in conjunction with methylation sites. The hybrid assembly bash script is designed to process Illumina paired-end short reads along with ONT long reads, and the long read assembly bash script only processes long reads. These two scripts will perform QC analysis on your reads, assemble the genomes, polish the genomes, perform QC on the genomes, and annotate the genomes. The methylation feature annotation python script is designed to process methylation output from MicrobeMod and genbank genome annotations from Prokka. The link to the MicrobeMod github is (https://github.com/cultivarium/MicrobeMod), the pipeline outlined was followed exactly to generate methylation output using these assemblies.

# Data
The Illumina short read paired end FASTQ files are available at:

The Oxford Nanopore Technologies long read FASTQ files are available at:

The Oxford Nanopore Technologies methylation aware unaligned BAM files are available at:

