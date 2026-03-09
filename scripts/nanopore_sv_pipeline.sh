#!/bin/bash
set -euo pipefail

# ======================================================
# Nanopore Long-Read Structural Variant Detection Pipeline
# Dataset: SRR15058164 (PRJNA744329)
# Author: Agrima Mateyal
# ======================================================

echo "======================================="
echo " Nanopore Structural Variant Pipeline "
echo "======================================="

# -------------------------------
# Define Project Directories
# -------------------------------

BASE_DIR=$(pwd)

RAW_DIR=$BASE_DIR/data/raw
FILTERED_DIR=$BASE_DIR/data/filtered
QC_DIR=$BASE_DIR/results/qc
ALIGN_DIR=$BASE_DIR/results/alignment
VARIANT_DIR=$BASE_DIR/results/variants
REF_DIR=$BASE_DIR/reference

mkdir -p $RAW_DIR $FILTERED_DIR $QC_DIR $ALIGN_DIR $VARIANT_DIR $REF_DIR

# -------------------------------
# Tool Versions
# -------------------------------

echo "[INFO] Tool versions:"
minimap2 --version
samtools --version
sniffles --version
bcftools --version

# ======================================================
# STEP 1: DATA DOWNLOAD
# ======================================================

echo "[STEP 1] Fetching run information..."

esearch -db sra -query PRJNA744329 | efetch -format runinfo > PRJNA744329_runs.csv

echo "[STEP 1] Downloading SRA dataset..."

prefetch SRR15058164

echo "[STEP 1] Converting SRA to FASTQ..."

fasterq-dump SRR15058164.sra --threads 8 -O $RAW_DIR --progress

echo "[STEP 1] Compressing FASTQ..."

gzip $RAW_DIR/SRR15058164.fastq

# ======================================================
# STEP 2: INITIAL READ STATISTICS
# ======================================================

echo "[STEP 2] Generating read statistics..."

seqkit stats -a $RAW_DIR/SRR15058164.fastq.gz > $QC_DIR/seqkit_stats.txt

# ======================================================
# STEP 3: QUALITY CONTROL
# ======================================================

echo "[STEP 3] Calculating N50..."

zcat $RAW_DIR/SRR15058164.fastq.gz | awk 'NR%4==2{print length($0)}' \
| sort -nr | awk '{sum+=$1; a[NR]=$1} END{for(i=1;i<=NR;i++){c+=a[i]; if(c>=sum/2){print "N50 =",a[i]; exit}}}' > $QC_DIR/n50.txt

echo "[STEP 3] Running NanoPlot on raw reads..."

NanoPlot \
--fastq $RAW_DIR/SRR15058164.fastq.gz \
-o $QC_DIR/nanoplot_raw \
--N50 --plots hex dot kde \
--loglength \
--threads 4

# ======================================================
# STEP 4: READ FILTERING
# ======================================================

echo "[STEP 4] Filtering reads (Q≥7, length≥1000)..."

zcat $RAW_DIR/SRR15058164.fastq.gz \
| NanoFilt -q 7 -l 1000 \
| gzip > $FILTERED_DIR/SRR15058164_filtered.fastq.gz

echo "[STEP 4] Running NanoPlot on filtered reads..."

NanoPlot \
--fastq $FILTERED_DIR/SRR15058164_filtered.fastq.gz \
-o $QC_DIR/nanoplot_filtered \
--loglength \
--N50 \
--summary

# ======================================================
# STEP 5: REFERENCE GENOME DOWNLOAD
# ======================================================

echo "[STEP 5] Downloading GRCh38 reference genome..."

wget -O $REF_DIR/GRCh38.fa.gz \
ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip -f $REF_DIR/GRCh38.fa.gz

# ======================================================
# STEP 6: ALIGNMENT
# ======================================================

echo "[STEP 6] Aligning reads using Minimap2..."

minimap2 \
-ax map-ont \
-t 8 \
$REF_DIR/GRCh38.fa \
$FILTERED_DIR/SRR15058164_filtered.fastq.gz \
> $ALIGN_DIR/aligned.sam

# ======================================================
# STEP 7: BAM PROCESSING
# ======================================================

echo "[STEP 7] Converting SAM to BAM..."

samtools view -Sb $ALIGN_DIR/aligned.sam > $ALIGN_DIR/aligned.bam

echo "[STEP 7] Sorting BAM..."

samtools sort \
$ALIGN_DIR/aligned.bam \
-o $ALIGN_DIR/aligned_sorted.bam

echo "[STEP 7] Marking duplicates..."

samtools markdup \
-r \
$ALIGN_DIR/aligned_sorted.bam \
$ALIGN_DIR/aligned_markdup.bam

echo "[STEP 7] Indexing BAM..."

samtools index $ALIGN_DIR/aligned_markdup.bam

echo "[STEP 7] Generating alignment statistics..."

samtools flagstat \
$ALIGN_DIR/aligned_markdup.bam \
> $QC_DIR/alignment_flagstat.txt

# ======================================================
# STEP 8: STRUCTURAL VARIANT CALLING
# ======================================================

echo "[STEP 8] Running Sniffles2 for SV detection..."

sniffles \
--input $ALIGN_DIR/aligned_markdup.bam \
--vcf $VARIANT_DIR/SRR15058164.sv.vcf

# ======================================================
# STEP 9: VARIANT STATISTICS
# ======================================================

echo "[STEP 9] Generating SV statistics..."

grep -v "^#" $VARIANT_DIR/SRR15058164.sv.vcf \
| wc -l \
> $VARIANT_DIR/sv_count.txt

grep -v "^#" $VARIANT_DIR/SRR15058164.sv.vcf \
| cut -f 5 \
| sort \
| uniq -c \
> $VARIANT_DIR/sv_types.txt

# ======================================================
# STEP 10: VARIANT FILTERING
# ======================================================

echo "[STEP 10] Filtering variants by size (SVLEN ≥ 50)..."

bcftools view \
-i 'ABS(SVLEN)>=50' \
$VARIANT_DIR/SRR15058164.sv.vcf \
-o $VARIANT_DIR/SRR15058164.filtered_svlen.vcf

echo "[STEP 10] Filtering variants by read support (SUPPORT ≥ 10)..."

bcftools view \
-i 'SUPPORT>=10' \
$VARIANT_DIR/SRR15058164.sv.vcf \
-o $VARIANT_DIR/SRR15058164.filtered_support.vcf

# ======================================================
# STEP 11: RESULT EXTRACTION
# ======================================================

echo "[STEP 11] Extracting example variants..."

grep -v "^#" $VARIANT_DIR/SRR15058164.sv.vcf \
| head -n 5 \
> $VARIANT_DIR/example_svs.txt

# ======================================================
# PIPELINE COMPLETE
# ======================================================

echo "======================================="
echo " Pipeline Completed Successfully "
echo " Results stored in: $VARIANT_DIR "
echo "======================================="
