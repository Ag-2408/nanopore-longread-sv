# Pipeline Workflow

This pipeline performs structural variant detection using Oxford Nanopore long-read sequencing data.

Steps:

1. Download sequencing data from NCBI SRA
2. Convert SRA to FASTQ format
3. Generate read statistics using SeqKit
4. Perform quality control using NanoPlot
5. Filter reads using NanoFilt
6. Align reads to GRCh38 using Minimap2
7. Process BAM files using SAMtools
8. Detect structural variants using Sniffles2
9. Filter variants using BCFtools
