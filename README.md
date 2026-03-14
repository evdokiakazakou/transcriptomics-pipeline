# Transcriptomics Pipeline

A pipeline for RNA-Seq data analysis covering read alignment, BAM file processing,
and differential expression analysis. Developed as part of the project for the
Proteomics and Functional Genomics course of the MSc program "Applied Bioinformatics".

---

## Table of Contents

1. [Overview](#1-overview)
2. [Part 1 — Read Alignment with HISAT2](#2-part-1--read-alignment-with-hisat2)
   - [Setup](#21-setup)
   - [Alignment](#22-alignment)
   - [SAM to BAM Conversion](#23-sam-to-bam-conversion)
3. [Part 2 — Differential Expression Analysis with metaseqR2](#3-part-2--differential-expression-analysis-with-metaseqr2)
   - [Dataset](#31-dataset)
   - [Targets File](#32-targets-file)
   - [Run metaseqR2](#33-run-metaseqr2)
   - [Alignment Metrics](#34-alignment-metrics)
4. [Acknowledgements](#4-acknowledgements)

---

## 1. Overview

This pipeline covers two parts:

| Part | Description | Tools |
|------|-------------|-------|
| Part 1 | Align RNA-Seq reads to a reference genome | HISAT2, SAMtools |
| Part 2 | Quality control and differential expression analysis | metaseqR2 (R) |

The dataset used in Part 2 comes from a study examining the transcriptome of prostatic
biopsies from patients with advanced hormone-naive prostate cancer, comparing pre- and
post-docetaxel chemotherapy treatment.

Reference: [Tzelepi et al.](https://pubmed.ncbi.nlm.nih.gov/) — Pre vs Post-docetaxel
prostatic biopsies transcriptome profiling.

---

## 2. Part 1 — Read Alignment with HISAT2

HISAT2 is a fast and sensitive splice-aware aligner for mapping RNA-Seq reads to a
reference genome. It performs spliced alignment by default, meaning it can handle reads
that span exon-exon junctions.

### 2.1 Setup

Create a working directory and download all required files:

```bash
mkdir align && cd align

# Download the FASTQ file
wget http://epigenomics.fleming.gr/~panos/appbio/human.fastq.gz

# Download and uncompress the hg19 genome index
wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
tar -xvf hg19_genome.tar.gz

# Download and unzip HISAT2
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
     -O hisat2-2.2.1-Linux_x86_64.zip
unzip hisat2-2.2.1-Linux_x86_64.zip

# Test HISAT2 installation
./hisat2-2.2.1/hisat2 --help

# Download and compile SAMtools (if not already available)
wget https://github.com/samtools/samtools/releases/download/1.23/samtools-1.23.tar.bz2
tar -xvf samtools-1.23.tar.bz2
cd samtools-1.23
./configure --without-curses --disable-bz2 --disable-lzma && make
cd ..

# Test SAMtools installation
./samtools-1.23/samtools --help
```

### 2.2 Alignment

HISAT2 performs spliced alignment by default. The reference genome index (hg19) must
be downloaded and placed in an accessible location before running the aligner.

```bash
./hisat2-2.2.1/hisat2 \
    -x hg19_genome/genome \
    -U human.fastq.gz \
    -S human_aligned.sam
```

- `-x` path to the genome index (without file extension)
- `-U` input FASTQ file (single-end; use `-1` and `-2` for paired-end)
- `-S` output SAM file

HISAT2 prints alignment statistics to the screen after completion, including the
total alignment rate. Note these statistics as they are required for the project report.

### 2.3 SAM to BAM Conversion

HISAT2 outputs a SAM file. Convert it to the compressed BAM format and index it:

```bash
# Convert SAM to BAM
./samtools-1.23/samtools view -bS human_aligned.sam > human_aligned.bam

# Sort the BAM file by genomic position
./samtools-1.23/samtools sort human_aligned.bam -o human_aligned_sorted.bam

# Index the sorted BAM file
./samtools-1.23/samtools index human_aligned_sorted.bam
```

- `-bS` converts SAM (-S) to binary BAM format (-b)
- Sorting is required before indexing

---

## 3. Part 2 — Differential Expression Analysis with metaseqR2

metaseqR2 is an R package for RNA-Seq quality control and differential expression
analysis. It supports multiple normalization and statistical testing methods and
produces an interactive HTML report.

### 3.1 Dataset

The dataset consists of paired prostatic biopsy RNA-Seq samples:

| Condition | Description |
|-----------|-------------|
| Pre | Biopsy before docetaxel chemotherapy treatment |
| Post | Biopsy after docetaxel chemotherapy treatment |

The comparison of interest is: **Post vs Pre-docetaxel**.

### 3.2 Targets File

The targets file is a tab-separated file that describes the samples and their metadata.
It is required by metaseqR2 to know which BAM files correspond to which condition.

Because the reads are paired-end and forward-stranded, the `strandedness` column must
be set to `forward`.

Create a file named `targets.txt` with the following structure:

```
samplename    filename              condition    paired    strandedness
pre_1         pre_sample1.bam       pre          TRUE      forward
pre_2         pre_sample2.bam       pre          TRUE      forward
pre_3         pre_sample3.bam       pre          TRUE      forward
post_1        post_sample1.bam      post         TRUE      forward
post_2        post_sample2.bam      post         TRUE      forward
post_3        post_sample3.bam      post         TRUE      forward
```

- `samplename` a short unique name for each sample
- `filename` the BAM file name
- `condition` the experimental group (pre or post)
- `paired` TRUE for paired-end reads
- `strandedness` forward, reverse, or no — based on the library preparation protocol

### 3.3 Run metaseqR2

Open R and run the following command. The normalization and statistical testing method
are both set to `deseq2` as required.

```r
library(metaseqR2)

metaseqr2(
    sampleList    = readTargets("targets.txt"),
    contrast      = c("post_vs_pre"),
    libsizeList   = NULL,
    normalization = "deseq2",
    statistics    = "deseq2",
    qcPlots       = c(
        "mds", "biodetection", "countsbio", "saturation",
        "correl", "boxplot", "meandiff", "meanvar",
        "volcano", "mastat"
    ),
    figFormat     = "png",
    exportWhere   = "./metaseqr2_output",
    reportTitle   = "Post vs Pre Docetaxel - Differential Expression Analysis",
    org           = "hg19",
    refdb         = "ensembl",
    report        = TRUE
)
```

**Parameter explanations:**

| Parameter | Description |
|-----------|-------------|
| `sampleList` | Reads the targets file to identify samples and conditions |
| `contrast` | The comparison to perform (Post vs Pre) |
| `normalization` | Normalization method — `deseq2` |
| `statistics` | Statistical testing method — `deseq2` |
| `qcPlots` | QC plots to include in the report |
| `exportWhere` | Directory where results will be saved |
| `org` | Reference genome organism — hg19 (human) |
| `refdb` | Annotation database — Ensembl |
| `report` | Generate an interactive HTML report |

The output includes an HTML report with all QC plots, a list of differentially
expressed genes, and normalized count tables.

### 3.4 Alignment Metrics

Generate alignment metrics using the metaseqR2 database and the same targets file:

```r
library(metaseqR2)

bf <- loadAnnotation(
    org   = "hg19",
    refdb = "ensembl",
    type  = "gene"
)

metaSeqAlignStats(
    targets   = readTargets("targets.txt"),
    bamRanges = bf
)
```

This produces per-sample alignment statistics based on the BAM files and the
annotated gene regions from the metaseqR2 database.

---

## 4. Acknowledgements

This pipeline is based on project material from the
**Proteomics and Functional Genomics** course,
MSc program "Applied Bioinformatics",
International Hellenic University (IHU), Thessaloniki.

Project designed by **Dr. Panagiotis Moulos**
(Fleming Institute, moulos@fleming.gr).
