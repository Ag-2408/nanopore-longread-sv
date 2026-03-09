#!/bin/bash
# ======================================================
# Long-Read Nanopore SV Pipeline Automation (SRR15058164)
# Based on HPRC Nanopore Data (PRJNA744329)
# ======================================================

# -------- DAY 0: Environment Setup --------
echo "[INFO] Creating and activating conda environment..."
conda create -y -n nsp python=3.10
conda activate nsp

# -------- DAY 1: Data Download --------
echo "[INFO] Fetching run info..."
esearch -db sra -query PRJNA744329 | efetch -format runinfo > PRJNA744329_runs.csv

echo "[INFO] Downloading SRA file..."
prefetch SRR15058164

echo "[INFO] Converting SRA to FASTQ..."
fasterq-dump /home/SDA/SRR15058164.sra --threads 8 -O /home/SDA/ --progress

echo "[INFO] Compressing FASTQ..."
gzip /home/SDA/SRR15058164.fastq
mv /home/SDA/SRR15058164.fastq.gz /home/SDA/raw/

echo "[INFO] Running SeqKit stats..."
seqkit stats -a /home/SDA/raw/SRR15058164.fastq.gz > /home/SDA/raw/seqkit_stats.txt

# -------- DAY 2: QC & Filtering --------
echo "[INFO] Calculating N50..."
zcat /home/SDA/raw/SRR15058164.fastq.gz | awk 'NR%4==2{print length($0)}' \
| sort -nr | awk '{sum+=$1; a[NR]=$1} END{for(i=1;i<=NR;i++){c+=a[i]; if(c>=sum/2){print "N50 =",a[i]; exit}}}' > /home/SDA/qc/n50.txt

echo "[INFO] Running NanoPlot (raw reads)..."
NanoPlot --fastq /home/SDA/raw/SRR15058164.fastq.gz -o /home/SDA/qc/nanoplot_out --N50 --plots hex dot kde --loglength --threads 4

echo "[INFO] Filtering reads (min Q7, length ≥1000)..."
zcat /home/SDA/raw/SRR15058164.fastq.gz | NanoFilt -q 7 -l 1000 | gzip > /home/SDA/filtered/SRR15058164_filtered.fastq.gz

echo "[INFO] Running NanoPlot (filtered reads)..."
NanoPlot --fastq /home/SDA/filtered/SRR15058164_filtered.fastq.gz -o /home/SDA/qc/nanoplot_filtered --loglength --N50 --summary

# -------- DAY 3: Alignment --------
echo "[INFO] Downloading GRCh38 reference genome..."
wget -O /home/SDA/align/GRCh38.fa.gz ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip /home/SDA/align/GRCh38.fa.gz

echo "[INFO] Running Minimap2 alignment..."
minimap2 -ax map-ont -t 8 /home/SDA/align/GRCh38.fa /home/SDA/filtered/SRR15058164_filtered.fastq.gz > /home/SDA/align/aligned.sam

echo "[INFO] Converting SAM to BAM..."
samtools view -Sb /home/SDA/align/aligned.sam > /home/SDA/align/aligned.bam

echo "[INFO] Sorting BAM..."
samtools sort /home/SDA/align/aligned.bam -o /home/SDA/align/aligned_sorted.bam

echo "[INFO] Marking duplicates..."
samtools markdup -r /home/SDA/align/aligned_sorted.bam /home/SDA/align/aligned_markdup.bam

echo "[INFO] Indexing BAM..."
samtools index /home/SDA/align/aligned_markdup.bam

echo "[INFO] Running flagstat..."
samtools flagstat /home/SDA/align/aligned_markdup.bam > /home/SDA/qc/align_flagstat.txt

# -------- DAY 4: Structural Variant Calling --------
echo "[INFO] Running Sniffles2..."
sniffles --input /home/SDA/align/aligned_markdup.bam --vcf /home/SDA/variants/SRR15058164.sv.vcf

echo "[INFO] Counting variants..."
grep -v "^#" /home/SDA/variants/SRR15058164.sv.vcf | wc -l > /home/SDA/variants/sv_count.txt
grep -v "^#" /home/SDA/variants/SRR15058164.sv.vcf | cut -f 5 | sort | uniq -c > /home/SDA/variants/sv_types.txt

# -------- DAY 5: Filtering & Annotation --------
echo "[INFO] Filtering variants (SVLEN ≥50)..."
bcftools view -i 'ABS(SVLEN) >= 50' /home/SDA/variants/SRR15058164.sv.vcf -o /home/SDA/variants/SRR15058164.filtered.vcf

echo "[INFO] Filtering variants (SUPPORT ≥10)..."
bcftools view -i 'SUPPORT>=10' /home/SDA/variants/SRR15058164.sv.vcf -o /home/SDA/variants/SRR15058164.filtered.vcf

echo "[INFO] Checking VEP annotation file..."
grep -v '^##' /home/SDA/variants/vep | sed -n '1,20p' > /home/SDA/variants/vep_head.txt

# -------- DAY 6: Visualization --------
echo "[INFO] Extracting example SVs..."
zgrep -v "^#" /home/SDA/variants/SRR15058164.sv.vcf.gz | head -n 5 > /home/SDA/variants/example_svs.txt

echo "[INFO] Pipeline completed successfully!"
