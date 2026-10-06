# MastersDissertation2025_HIVDataAnalysis
This repository contains the python code used for my MSc Genomic Medicine dissertation (University of Oxford 2025) on *Quantitative analysis of suspected HIV proviral DNA reads observed from strand-specific HIV metagenomic sequencing data.*

**Summary**
HIV proviral DNA persists in cells even when treatment suppresses viral load, so it could help detect drug resistance early. This project tests whether Strand-Specific Metagenomic sequencing built for HIV RNA, the veSEQ-HIV method developed by the Bonsall group also captures proviral DNA [Bonsall et al. (2020)] (https://pubmed.ncbi.nlm.nih.gov/32669382/).

**The idea:** HIV RNA appears almost entirely as forward-strand reads, while double-stranded proviral DNA also produces reverse-strand reads. The first 454 bases of the proviral genome have no RNA equivalent, so reads there should split roughly 1:1 between strands if proviral DNA is present.

The pipeline in this repository separates reads by strand, computes per-base coverage and forward/reverse ratios, groups samples by viral load, and analyses reads that span the proviral-specific regions. 

This dissertation was conceptualized, written and submitted between May 2025 to July 2025

**Data and credits**
No sequencing data is included. All data belongs to the Bonsall group and was collected under the PANGEA-HIV consortium. Ethical approval restricted analysis to HIV sequences, so all non-HIV reads were removed before my analysis.
The starting point for every script below is a set of sorted and indexed BAM files (sorted and indexed with samtools).
Pre-processing used the Located K-mer Assembler (LKA), developed by Dr Nicholas Grayson. This assembler tool is not part of this repository and is available from the Bonsall group on request.

All scripts here were written by me, with supervision from Dr David Bonsall and Dr Nicholas Grayson.

## Workflow

### Workflow at a glance

| Step | What it does | Script | Input | Main output |
|------|--------------|--------|-------|-------------|
| 1 | Strand-specific coverage, pooled PacBio data | `LogCoverage_and_Coverage_ratios.py` | Folder of sorted, indexed BAM files | Coverage and ratio plots, `strand_summary.txt`, `coverage_table.tsv` |
| 2a | Percentage of forward vs reverse reads per pool | `Summary_of_Pacbio_read_percentages.py` | `strand_summary.txt` files from Step 1 | Stacked bar plot |
| 2b | Forward vs reverse read counts per sample | `Forward_to_reverse_read_strand_per_sample.py` | Per-sample `strand_summary.txt` files | `strand_log_scatter.png` |
| 3 | Strand-specific coverage, individual Illumina data | `PerSample_LogCoverage_and_Coverage_ratios.py` | Folder of sorted, indexed BAM files | Per-sample plots and tables, summary TSV, interactive HTML plot |
| 4 | Reads spanning a chosen HXB2 region | `Analysis_of_reads_spanning_fiveprime_U3_R_region.py` | One or more sorted BAM files | Summary CSV, plots, optional reverse-read files |

### Workflow steps listed

The starting data for my analysis consisted of BAM files. The steps of my workflow that involved custom Python scripts are as follows:

1. Strand-specific BAM Coverage Analysis of Pacbio datasets (Code: LogCoverage_and_Coverage_ratios.py)

2.  (a)Summary plot of read percentages from pooled datasets (Code: Summary_of_Pacbio_read_percentages.py)

   (b)Log-scaled scatter plot of Forward vs Reverse read counts per sample fo individual datasets (Code: Forward_to_reverse_read_strand_per_sample.py)
   
3. Strand-specific BAM Coverage Analysis of individual Illumina datasets (Code: PerSample_LogCoverage_and_Coverage_ratios.py)
   
4. Strand-specific analysis of reads that span user defined HXB2 genomic regions (Code: Analysis_of_reads_spanning_fiveprime_U3_R_region.py)


The codes I used for each step of my analysis are documented here for reproducibility
   

# Step 1: Strand-Specific BAM Coverage Analysis of pooled datasets

Python pipeline for separating BAM reads by strand, computing per-base coverage, and visualizing strand bias in HIV sequencing datasets (for example, HIV proviral integration studies).
See LogCoverage_and_Coverage_ratios.py for code used

## Key Features

- Separates forward and reverse reads into distinct BAM files.
- Computes per-base genome-wide coverage (log10 scaled).
- Plots strand-specific coverage profiles.
- Calculates and visualizes the forward/reverse strand coverage ratio.
- Outputs tabulated per-position coverage and ratios.
- Fully automated for batch processing multiple samples.

## Why Use This?

In strand-specific studies (e.g. HIV proviral detection, transcriptomics), it's critical to distinguish the orientation of aligned reads. 
This analysis enables quick, reproducible inspection of strand bias across reference genomes.

---

## Requirements
- Python 3.7+
- `pysam`, `numpy`, `matplotlib.pyplot`

Install requirements via pip:

```bash
pip install pysam matplotlib numpy
```

## Input

- A directory of **sorted and indexed BAM files**

## Output
- A log coverage plot showing coverage of forward reads (marked blue) versus reverse reads (marked red)
- A coverage ratio graph fowing the forward:reverse strand ratios at each position across the genome
- A `strand_summary.txt` file and a `coverage_table.tsv` showing the ratio of forward:reverse reads at each position in the genome

## Usage

```bash
python LogCoverage_and_Coverage_ratios.py --input_dir bam_folder --output_dir coverage_outputs


Install requirements via pip:

pip install matplotlib numpy
```
# Step 2(a): Summary bar plot of read percentages from pooled Pacbio datasets

This script visualizes the proportion of **forward vs. reverse strand reads** across different sample sets based on strand-specific summary statistics.

It parses `strand_summary.txt` files (generated by a prior BAM processing script) and produces a **stacked bar plot** showing the percentage of reads aligned to each strand per sample.

See Summary_of_Pacbio_read_percentages.py for the code used

---
## Requirements

- Python 3.7 and above
- `matplotlib`
- `numpy`

Install requirements via pip:

```bash
pip install matplotlib numpy
```
## Input

-  `strand_summary.txt` files

  
## Output

Each sample is represented by a bar divided into:
- **Forward strand reads** (top, blue)
- **Reverse strand reads** (bottom, red)
Percentages are displayed on the bars for easy interpretation.
---

# Step 2(b): Log-scaled scatter plot of Forward vs Reverse read counts per sample fo individual datasets

This Python script generates a **log-log scatter plot** of **forward vs reverse read counts** extracted from `strand_summary.txt` files produced during PacBio sequencing analysis.

It is ideal for identifying strand bias, visualizing sample coverage, and spotting potential sequencing anomalies across many individual samples.

---
## Features

- Parses `strand_summary.txt` files from nested subdirectories  
- Extracts forward and reverse strand read counts  
- Generates a **log–log scatter plot**  
- Includes a diagonal **reference line (F = R)**  
- Clean, minimalist layout with:
- No grid clutter
- Legend positioned at the **top left**

## Expected folder structure
```
coverage_output/
  Sample_01/
      forward_reads.bam
      forward_reads.sorted.bam
      forward_reads.sorted.bam.bai
      reverse_reads.bam
      reverse_reads.sorted.bam
      reverse_reads.sorted.bam.bai
      strand_summary.txt
   Sample_02/
      forward_reads.bam
      forward_reads.sorted.bam
      forward_reads.sorted.bam.bai
      reverse_reads.bam
      reverse_reads.sorted.bam
      reverse_reads.sorted.bam.bai
      strand_summary.txt

Each strand summary.txt should contain lines like:
↳ Forward strand: 29124 (46.84%)
↳ Reverse strand: 33054 (53.16%)

```   

## Usage

```bash
python3 Forward_to_reverse_strand_per_sample.py --summary_root .. --output strand_log_scatter.png

```
## Output

A scatter plot named strand_log_scatter.png (or your specified --output) showing:
   X-axis: forward read counts (log scale)
   Y-axis: reverse read counts (log scale)
   Diagonal y = x reference line
   Legend in the top-left corner

## Requirements

- Python 3.7 and above
- `matplotlib`
- `numpy`
- `re`

Install requirements via pip:

```bash
pip install matplotlib numpy
```
   
# Step 3: Strand-Specific BAM Coverage Analysis of individual datasets

Python pipeline that processes BAM files to: 

- Separate reads by strand (forward vs reverse)
- Compute and visualize strand-specific coverage
- Calculates and visualize the forward/reverse strand coverage ratio to assess strand balance
- Generate log-scaled interactive and static plots for forward vs reverse strand coverage

See PerSample_LogCoverage_and_Coverage_ratios.py for code used

 ## Features
- Strand separation using `pysam`
- Optional strand flipping for dUTP libraries
- Per-base coverage calculation
- Log-scaled coverage plots with dropout regions marked
- Forward/Reverse ratio analysis with strand ratio plots
- Summary statistic tsv files
- Interactive HTML visualization plots via plotly

 ## Requirements

- Python 3.7+
- `pysam`, `numpy`, `matplotlib`, `seaborn`, `pandas`, `plotly`

Install dependencies:
```bash
pip install pysam numpy matplotlib seaborn pandas plotly
```

## Usage

```bash
python PerSample_LogCoverage_and_Coverage_ratios.py \
    --input_dir ./bam_files \
    --output_dir ./results \
    --flip_strands #optional for dUTP strand flipped libraries
```

## Input

- A directory of **sorted and indexed BAM files**

## Output
- A log coverage plot showing coverage of forward reads (marked blue) versus reverse reads (marked red)
- A coverage ratio graph fowing the forward:reverse strand ratios at each position across the genome
- A `strand_summary.txt` file and a `coverage_table.tsv` showing the ratio of forward:reverse reads at each position in the genome
  
 ### Example output
```
coverage_output/
  Sample_01/
      forward_reads.bam
      forward_reads.sorted.bam
      forward_reads.sorted.bam.bai
      reverse_reads.bam
      reverse_reads.sorted.bam
      reverse_reads.sorted.bam.bai
      sample_01_combined_log10_coverage.png
      sample_01_coverage_ratio.png
      sample_01_coverage_table.tsv
   Sample_02/

strand_summary.tsv          
strand_scatter_log.png      
interactive_scatter.html
```   




# Step 4: Strand-specific analysis of reads that span user defined HXB2 genomic regions

A Python pipeline that analyzes strand-specific reads that span a user-defined genomic region from one or more BAM files. It supports visual summaries, reverse strand extraction, and publication-ready plots.

---

## Features

- Counts reads spanning a defined region (default: 450–633 bp).
- Classifies reads as forward or reverse strand based on alignment orientation.
- Outputs:
  - Read summary CSV file
  - Percentage Forward/Reverse count bar plots
  - Grouped bar plots (optional)
  - Lollipop plots for % reverse strand reads (optional)
  - Reverse strand FASTA, BAM, and TSV output (if `--reverse-only` is set)

---

## Requirements

- Python 3.6+
- `pysam`  `matplotlib`  `pandas` 
- `argsparse`  `csv`  `os`  `re` (standard libraries)

Install required libraries:
```bash
pip install matplotlib pandas
```
## Usage

```bash

python Analysis_of_reads_spanning_fiveprime_U3_R_region.py *.bam -o output_folder --region 450 633 --plot-grouped --plot-lollipop

```
### Required arguments:
*.bam: One or more sorted BAM files.

### Optional arguments:
```
-o, --output_dir: Output directory (default: strand_read_output)

--region START END: Genomic region to analyze (default: 450–633)

--reverse-only: Only process and export reverse strand reads

--plot-grouped: Include grouped bar chart of forward/reverse counts

--plot-lollipop: Include lollipop plot for % reverse reads

```
## Input

- A directory of **sorted and indexed BAM files**

## Example Output
```
strand_read_output/
   sample1_spanning_reads.tsv
   sample1_reverse_reads.fasta
   sample1_reverse_reads.bam
   strand_read_summary.csv
   strand_read_summary.png
   grouped_read_counts.png
   lollipop_reverse_percent.png
```

