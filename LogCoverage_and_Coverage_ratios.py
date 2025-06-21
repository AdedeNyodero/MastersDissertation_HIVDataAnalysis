#!/usr/bin/env python3

import pysam
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys
import glob


def separate_strands(input_bam, forward_bam_path, reverse_bam_path, summary_path):
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    forward_out = pysam.AlignmentFile(forward_bam_path, "wb", header=bamfile.header)
    reverse_out = pysam.AlignmentFile(reverse_bam_path, "wb", header=bamfile.header)

    forward_count = reverse_count = unmapped_count = total_count = 0

    for read in bamfile:
        total_count += 1
        if read.is_unmapped:
            unmapped_count += 1
            continue
        if read.is_reverse:
            reverse_out.write(read)
            reverse_count += 1
        else:
            forward_out.write(read)
            forward_count += 1

    bamfile.close()
    forward_out.close()
    reverse_out.close()

    mapped_count = forward_count + reverse_count
    summary = (
        "--Strand Read Summary--\n"
        f"Total reads:          {total_count}\n"
        f"Unmapped reads:       {unmapped_count}\n"
        f"Mapped reads:         {mapped_count}\n"
        f"  ↳ Forward strand:   {forward_count} ({forward_count/mapped_count*100:.2f}%)\n"
        f"  ↳ Reverse strand:   {reverse_count} ({reverse_count/mapped_count*100:.2f}%)\n"
    )

    print(f"\n[Sample: {os.path.basename(input_bam)}] {summary}")
    with open(summary_path, "w") as f:
        f.write(summary)


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

    positions_to_mark = {
        454: 'black',
        634: 'grey',
        9086: 'grey',
        9540: 'black',
        9649: 'black'
    }

    for ref in refs:
        cov_f = forward_cov[ref]
        cov_r = reverse_cov[ref]

        log_f = np.log10(cov_f + 1)
        log_r = np.log10(cov_r + 1)

        plt.figure(figsize=(12, 4))
        plt.plot(log_f, color='blue', label='Forward strand')
        plt.plot(log_r, color='red', label='Reverse strand')

        yticks = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
        plt.yticks(yticks, [f"{t:.1f}" for t in yticks])
        plt.ylim(1.0, 7.0)
        plt.ylabel("Log10-scaled coverage")

        plt.xlabel("Genomic Position")
        plt.title(f"Strand Coverage (log10 scale) for {sample_name}")

        length = len(cov_f)
        xticks = np.arange(0, length + 1, 500)
        plt.xticks(xticks)

        _, y_max = plt.ylim()
        for xpos, color in positions_to_mark.items():
            plt.axvline(x=xpos, color=color, linestyle=':', linewidth=1)
            plt.text(xpos, y_max * 0.95, f"Pos {xpos}", rotation=90,
                     va='top', ha='center', fontsize=8, color=color, clip_on=True)


        plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
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

    positions_to_mark = {
        454: 'black',
        634: 'grey',
        9086: 'grey',
        9540: 'black',
        9649: 'black'
    }

    for ref in refs:
        arr = ratio_dict[ref]
        plt.figure(figsize=(12, 4))
        plt.plot(arr, color='purple', label='Forward/Reverse ratio')
        plt.xlabel('Genomic Position')
        plt.ylabel('Coverage ratio (F/R)')
        plt.title(f'Coverage Ratio for {sample_name}')

        length = len(arr)
        xticks = np.arange(0, length + 1, 500)
        plt.xticks(xticks)

        _, y_max = plt.ylim()
        for xpos, color in positions_to_mark.items():
            plt.axvline(x=xpos, color=color, linestyle=':', linewidth=1)
            plt.text(xpos, y_max * 0.95, f"{xpos}", rotation=90,
                     va='top', ha='center', fontsize=8, color=color, clip_on=True)

        # Legend outside the plot area
        plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
        plt.tight_layout()
        out_png = os.path.join(output_dir, f"{ref}_coverage_ratio.png")
        plt.savefig(out_png, bbox_inches='tight')
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
            out.write("Position\t(f+1)\t(r+1)\tRatio\n")
            for pos in positions:
                out.write(f"{pos}\t{f[pos]+1}\t{r[pos]+1}\t{ratios[pos]:.2f}\n")
        print(f"Saved coverage table: {tsv_path}")


def sort_and_index_bam(in_bam, out_bam):
    pysam.sort("-o", out_bam, in_bam)
    pysam.index(out_bam)


def main():
    parser = argparse.ArgumentParser(
        description="Batch process all BAMs in a folder and plot coverage, ratio, and save tables."
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing sorted .bam files")
    parser.add_argument("--output_dir", default="coverage",
                        help="Directory for all coverage outputs")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    bam_files = sorted(glob.glob(os.path.join(args.input_dir, "*.bam")))
    if not bam_files:
        print(f"No BAM files found in {args.input_dir}")
        sys.exit(1)

    for bam in bam_files:
        sample = os.path.basename(bam).rsplit('.bam', 1)[0].split('_B')[0]
        sample_out = os.path.join(args.output_dir, sample)
        os.makedirs(sample_out, exist_ok=True)

        f_bam = os.path.join(sample_out, "forward_reads.bam")
        r_bam = os.path.join(sample_out, "reverse_reads.bam")
        summary = os.path.join(sample_out, "strand_summary.txt")

        print(f"\n=== Processing {sample} ===")
        separate_strands(bam, f_bam, r_bam, summary)

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

    print("\n All samples processed.")


if __name__ == "__main__":
    main()