# masters_dissertation_data_analysis_25
This repository contains code used for my masters dissertation on Quantitative analysis of suspected HIV proviral DNA reads observed from the strand-specific HIV RNA metagenomic libraries. 
This dissertation was conceptualized, writtedn and submitted between May 2025 to July 2025

The steps of my work that used custon Python code are:
1. Strand-specific BAM Coverage Analysis (Code: LogCoverage_and_Coverage_ratios.py)





The python codes I used for each step of my analysis is documented here for reproducibility
   

# Step 1: Strand-Specific BAM Coverage Analysis

Python pipeline for separating BAM reads by strand, computing per-base coverage, and visualizing strand bias in HIV sequencing datasets (for example, HIV proviral integration studies).

  # Key Features

- Separates forward and reverse reads into distinct BAM files.
- Computes per-base genome-wide coverage (log₁₀-scaled).
- Plots strand-specific coverage profiles.
- Calculates and visualizes the forward/reverse strand coverage ratio.
- Outputs tabulated per-position coverage and ratios.
- Fully automated for batch processing multiple samples.

  # Why Use This?

In strand-specific studies (e.g. HIV proviral detection, transcriptomics), it's critical to distinguish the orientation of aligned reads. 
This tool enables quick, reproducible inspection of strand bias across reference genomes.

---

  # Requirements

- Python version 3.12.3

  # Dependencies
- pysam (https://pysam.readthedocs.io/)
- matplotlib.pyplot
- numpy

  Install dependencies with:
pip install pysam matplotlib numpy (on bash)

  # Input

- A directory of **sorted and indexed BAM files**

  # Usage

```bash
python LogCoverage_and_Coverage_ratios.py --input_dir bam_folder --output_dir coverage_outputs

 
