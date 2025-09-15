# Bacterial-Methylation
This repository contains two bash scripts for bacterial genome assembly, one bash script for methylome assembly, and one lsf file for methylaiton calling using Dorado. The hybrid assembly bash script is designed to process Illumina paired-end short reads along with ONT long reads, and the long read assembly bash script only processes long reads. These two scripts will perform QC analysis on your reads, assemble the genomes, polish the genomes, perform QC on the genomes, and annotate the genomes. We used the tool MicrobeMod to get methylation site calls, motifs, and restriction modification gene output from our bacterial strains. The link to the MicrobeMod github is (https://github.com/cultivarium/MicrobeMod) and we followed the pipeline exactly. The methylation calling lsf file is designed to process raw ONT pod5 files into demultiplexed, unaligned .bam files and then rename them from the barcode number to the replicate ID.

# Data
The Illumina short read PE FastQ files, ONT long read FastQ files, and ONT uBAM methylation files are available at: https://doi.org/10.5281/zenodo.15555625 

The raw ONT pod5 files are available at: 

The preprint that this anaysis was performed in is at: doi: https://doi.org/10.1101/2025.06.04.657942 
