# Nanopore Long-Read Structural Variant Pipeline

![Bioinformatics](https://img.shields.io/badge/domain-bioinformatics-blue)
![Pipeline](https://img.shields.io/badge/pipeline-long--read%20SV%20analysis-green)
![Language](https://img.shields.io/badge/language-bash-yellow)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

A bioinformatics workflow for detecting **structural variants (SVs)** from **Oxford Nanopore long-read sequencing data**.

This pipeline automates processing from **raw sequencing data → structural variant discovery → filtered results**.

Dataset used: **SRR15058164 (PRJNA744329)**

---

# Overview

Long-read sequencing technologies such as **Oxford Nanopore** allow accurate detection of **large structural variants** that are difficult to identify with short-read sequencing.

This pipeline performs:

* Data download from NCBI SRA
* FASTQ generation
* Read statistics and quality control
* Read filtering
* Alignment to the human reference genome
* Structural variant detection
* Variant filtering and summary generation

---

# Pipeline Workflow

```mermaid
flowchart TD

A[SRA Dataset]

B[FASTQ Generation<br>SRA Toolkit]

C[Read Statistics<br>SeqKit]

D[Quality Control<br>NanoPlot]

E[Read Filtering<br>NanoFilt]

F[Alignment to GRCh38<br>Minimap2]

G[BAM Processing<br>SAMtools]

H[Structural Variant Detection<br>Sniffles2]

I[Variant Filtering<br>BCFtools]

J[Structural Variant Results]


A --> B
B --> C
C --> D
D --> E
E --> F
F --> G
G --> H
H --> I
I --> J


classDef data fill:#1f77b4,color:white
classDef qc fill:#ff7f0e,color:white
classDef align fill:#9467bd,color:white
classDef sv fill:#2ca02c,color:white


class A data
class B,C,D,E qc
class F,G align
class H,I,J sv
```
---

# Tools Used

| Tool        | Purpose                      |
| ----------- | ---------------------------- |
| SRA Toolkit | Download sequencing data     |
| SeqKit      | FASTQ statistics             |
| NanoPlot    | Read quality analysis        |
| NanoFilt    | Filtering low-quality reads  |
| Minimap2    | Long-read alignment          |
| SAMtools    | BAM processing               |
| Sniffles2   | Structural variant detection |
| BCFtools    | Variant filtering            |

---

# Repository Structure

```
nanopore-longread-sv
│
├── scripts
│   └── nanopore_sv_pipeline.sh
│
├── environment
│   └── environment.yml
│
├── data
│
├── results
│
└── README.md
```

---

# Installation

Clone the repository:

```
git clone https://github.com/Ag-2408/nanopore-longread-sv.git
cd nanopore-longread-sv
```

Create the conda environment:

```
conda env create -f environment/environment.yml
conda activate nanopore_sv_pipeline
```

---

# Running the Pipeline

Execute the pipeline script:

```
bash scripts/nanopore_sv_pipeline.sh
```

---

# Output

The pipeline generates:

* FASTQ statistics
* Quality control reports
* Alignment statistics
* Structural variant VCF files
* Filtered structural variants
* Variant summary reports

All results are stored in the **results/** directory.

---

# Author

**Agrima Mateyal**
MSc Bioinformatics
