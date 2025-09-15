# Bacterial-Methylation
This repository contains two bash scripts for bacterial genome assembly, one lsf file for methylation calling using ONT's Dorado basecaller, and one bash script for methylome assembly. 

Genome Assembly: The hybrid assembly bash script is designed to process Illumina paired-end short reads along with ONT long reads, and the long read assembly bash script only processes long reads. These two scripts will perform QC analysis on your reads, assemble the genomes, polish the genomes, perform QC on the genomes, and annotate the genomes. 

Methylation Calling: The methylation calling lsf file is designed to process raw ONT pod5 files into demultiplexed, unaligned .bam files and then rename them from the barcode number to the replicate ID. We used the tool MicrobeMod to get methylation site calls, motifs, and restriction modification gene output from our bacterial strains. The link to the MicrobeMod github is (https://github.com/cultivarium/MicrobeMod). The bash script for methylome assembly is designed to take both hybrid and long-read only assemblies in .fasta format and run the MicrobeMod pipeline on them. This includes alignment of methylation uBAM files to the reference hybrid or long-read only genome, running the MicrobeMod call methylation pipeline, and running the MicrobeMod annotate restriction modification systems pipeline.

# Data Availibility
The Illumina short read PE FastQ files, ONT long read FastQ files, and ONT uBAM methylation files are available at: https://doi.org/10.5281/zenodo.15555625 

The raw ONT pod5 files are available at: 

The preprint that this anaysis was performed in is at: doi: https://doi.org/10.1101/2025.06.04.657942 
