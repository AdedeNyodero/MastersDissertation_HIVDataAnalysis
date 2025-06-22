#!/usr/bin/env python3

import pysam
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys
import glob
import pandas as pd
import seaborn as sns
import plotly.express as px

def separate_strands(input_bam, forward_bam_path, reverse_bam_path, flip_strands=False):
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    forward_out = pysam.AlignmentFile(forward_bam_path, "wb", header=bamfile.header)
    reverse_out = pysam.AlignmentFile(reverse_bam_path, "wb", header=bamfile.header)

    forward_count = reverse_count = unmapped_count = total_count = 0

    for read in bamfile:
        total_count += 1
        if read.is_unmapped:
            unmapped_count += 1
            continue

        is_forward = not read.is_reverse
        if flip_strands:
            is_forward = not is_forward

        if is_forward:
            forward_out.write(read)
            forward_count += 1
        else:
            reverse_out.write(read)
            reverse_count += 1

    bamfile.close()
    forward_out.close()
    reverse_out.close()

    return forward_count, reverse_count, unmapped_count, total_count

def check_index_exists(bam_path):
    bai = bam_path + ".bai"
    if not os.path.exists(bai):
        raise FileNotFoundError(f"Missing BAM index: {bai}. Ensure BAM is sorted & indexed.")

def compute_coverage(bam_path):
    check_index_exists(bam_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    coverage = {}
    for ref in bam.references:
        length = bam.get_reference_length(ref)
        cov = np.zeros(length, dtype=int)
        for col in bam.pileup(ref, stepper="all", truncate=True, max_depth=1_000_000_000):
            pos = col.reference_pos
            if pos < length:
                cov[pos] = col.nsegments
        coverage[ref] = cov
    bam.close()
    return coverage

def plot_combined_coverage_log(forward_cov, reverse_cov, output_dir, sample_name, refs=None):
    if refs is None:
        refs = sorted(forward_cov.keys())

    for ref in refs:
        cov_f = forward_cov[ref]
        cov_r = reverse_cov[ref]

        log_f = np.log10(cov_f + 1)
        log_r = np.log10(cov_r + 1)

        plt.figure(figsize=(14, 5))
        plt.plot(log_f, color='blue', label='Forward strand')
        plt.plot(log_r, color='red', label='Reverse strand')

        # Overlay low reverse coverage dropouts as black dots
        dropout_positions = np.where(log_r < 1.2)[0]
        plt.scatter(dropout_positions, log_r[dropout_positions], color='black', s=5, alpha=0.6, label='Reverse < 1.2')

        plt.yticks(np.arange(0, 7.5, 0.5))
        plt.ylim(0.5, 7)
        plt.ylabel("Log10-scaled coverage")
        plt.xlabel("Genomic Position")
        plt.title(f"Strand Coverage (log10 scale) for {sample_name}")
        xticks = np.arange(0, len(log_f) + 1, 500)
        plt.xticks(xticks)
        plt.legend()
        plt.tight_layout()
        out_png = os.path.join(output_dir, f"{ref}_combined_log10_coverage.png")
        plt.savefig(out_png, bbox_inches='tight')
        plt.close()

def compute_strand_ratio(forward_cov, reverse_cov):
    ratio = {}
    for ref in forward_cov:
        f = forward_cov[ref]
        r = reverse_cov[ref]
        ratio[ref] = (f + 1) / (r + 1)
    return ratio

def plot_ratio(ratio_dict, output_dir, sample_name, refs=None):
    if refs is None:
        refs = sorted(ratio_dict.keys())

    for ref in refs:
        arr = ratio_dict[ref]
        plt.figure(figsize=(12, 4))
        plt.plot(arr, color='purple', label='Forward/Reverse ratio')
        plt.xlabel('Genomic Position')
        plt.ylabel('Coverage ratio (F/R)')
        plt.title(f'Coverage Ratio for {sample_name}')
        plt.legend(loc='upper right')
        plt.tight_layout()
        out_png = os.path.join(output_dir, f"{ref}_coverage_ratio.png")
        plt.savefig(out_png)
        plt.close()

def save_coverage_table(forward_cov, reverse_cov, ratio_dict, output_dir, sample_name, refs=None):
    if refs is None:
        refs = sorted(forward_cov.keys())
    for ref in refs:
        f = forward_cov[ref]
        r = reverse_cov[ref]
        ratios = ratio_dict[ref]
        positions = np.arange(len(f))
        tsv_path = os.path.join(output_dir, f"{ref}_coverage_table.tsv")
        with open(tsv_path, 'w') as out:
            out.write("Position	(f+1)	(r+1)	Ratio")
            for pos in positions:
                out.write(f"{pos}	{f[pos]+1}	{r[pos]+1}	{ratios[pos]:.2f}")

def sort_and_index_bam(in_bam, out_bam):
    pysam.sort("-o", out_bam, in_bam)
    pysam.index(out_bam)

def main():
    parser = argparse.ArgumentParser(
        description="Master pipeline for strand separation, coverage plotting, and interactive QC."
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing sorted .bam files")
    parser.add_argument("--output_dir", default="coverage_output", help="Directory for all outputs")
    parser.add_argument("--flip_strands", action="store_true", help="Flip strand orientation (e.g., dUTP libraries)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    bam_files = sorted(glob.glob(os.path.join(args.input_dir, "*.bam")))
    if not bam_files:
        print(f"No BAM files found in {args.input_dir}")
        sys.exit(1)

    summary_data = []

    for bam in bam_files:
        sample = os.path.basename(bam).rsplit('.bam', 1)[0].split('_B')[0]
        sample_out = os.path.join(args.output_dir, sample)
        os.makedirs(sample_out, exist_ok=True)

        f_bam = os.path.join(sample_out, "forward_reads.bam")
        r_bam = os.path.join(sample_out, "reverse_reads.bam")

        print(f"\n=== Processing {sample} ===")
        forward_count, reverse_count, unmapped_count, total_count = separate_strands(
            bam, f_bam, r_bam, flip_strands=args.flip_strands
        )

        sf = f_bam.replace(".bam", ".sorted.bam")
        sr = r_bam.replace(".bam", ".sorted.bam")
        sort_and_index_bam(f_bam, sf)
        sort_and_index_bam(r_bam, sr)

        cov_f = compute_coverage(sf)
        cov_r = compute_coverage(sr)

        plot_combined_coverage_log(cov_f, cov_r, sample_out, sample)
        ratios = compute_strand_ratio(cov_f, cov_r)
        plot_ratio(ratios, sample_out, sample)
        save_coverage_table(cov_f, cov_r, ratios, sample_out, sample)

        summary_data.append({
            "Sample": sample,
            "Forward_Reads": forward_count,
            "Reverse_Reads": reverse_count,
            "Total_Reads": total_count,
            "Unmapped_Reads": unmapped_count,
            "Mapped_Reads": forward_count + reverse_count,
            "Percent_Forward": 100 * forward_count / (forward_count + reverse_count) if (forward_count + reverse_count) > 0 else 0,
            "Log2_FR": np.log2(forward_count / reverse_count) if forward_count > 0 and reverse_count > 0 else np.nan,
            "Flip_Strand_Applied": args.flip_strands
        })

    # Save CSV summary
    df_summary = pd.DataFrame(summary_data)
    csv_path = os.path.join(args.output_dir, "strand_summary.csv")
    df_summary.to_csv(csv_path, index=False)

    # Static log-scale scatter plot
    plt.figure(figsize=(8, 8))
    sns.scatterplot(data=df_summary, x="Forward_Reads", y="Reverse_Reads")
    plt.xscale("log")
    plt.yscale("log")
    lims = [max(df_summary[["Forward_Reads", "Reverse_Reads"]].min().min(), 1),
            df_summary[["Forward_Reads", "Reverse_Reads"]].max().max()]
    plt.plot(lims, lims, color='gray', linestyle='--', label='F=R line')
    plt.xlabel("Forward Read Count (log scale)")
    plt.ylabel("Reverse Read Count (log scale)")
    plt.title("Log-Scaled Forward vs Reverse Reads per Sample")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "strand_scatter_log.png"))

    # Interactive Plotly output
    fig = px.scatter(
        df_summary,
        x="Forward_Reads",
        y="Reverse_Reads",
        hover_name="Sample",
        log_x=True,
        log_y=True,
        title="Interactive Log-Scaled Scatter Plot of Forward vs Reverse Reads",
        labels={"Forward_Reads": "Forward Read Count", "Reverse_Reads": "Reverse Read Count"}
    )
    fig.add_shape(
        type='line',
        x0=1, y0=1,
        x1=df_summary["Forward_Reads"].max(),
        y1=df_summary["Forward_Reads"].max(),
        line=dict(color='gray', dash='dash')
    )
    fig.write_html(os.path.join(args.output_dir, "interactive_scatter.html"))

    print(f"\n All samples processed. Results in: {args.output_dir}")

if __name__ == "__main__":
    main()
