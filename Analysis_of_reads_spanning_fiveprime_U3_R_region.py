import pysam
import os
import argparse
import matplotlib.pyplot as plt
import csv
import re
import pandas as pd

def analyze_bam(bam_path, output_dir, region_start, region_end, reverse_only=False):
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    tsv_path = os.path.join(output_dir, f"{sample_name}_spanning_reads.tsv")
    reverse_tsv_path = os.path.join(output_dir, f"{sample_name}_reverse_reads.tsv")
    reverse_fasta_path = os.path.join(output_dir, f"{sample_name}_reverse_reads.fasta")
    reverse_bam_path = os.path.join(output_dir, f"{sample_name}_reverse_reads.bam")

    bamfile = pysam.AlignmentFile(bam_path, "rb")
    header = bamfile.header
    reverse_bam_out = pysam.AlignmentFile(reverse_bam_path, "wb", header=header) if reverse_only else None
    fasta_out = open(reverse_fasta_path, "w") if reverse_only else None
    tsv_out = open(reverse_tsv_path, "w") if reverse_only else open(tsv_path, "w")
    tsv_out.write("QueryName\tStrand\tRefStart\tRefEnd\n")

    total = forward = reverse = 0
    for read in bamfile.fetch():
        if read.is_unmapped:
            continue
        if read.reference_start < region_start and read.reference_end > region_end:
            strand = "reverse" if read.is_reverse else "forward"
            if reverse_only:
                if read.is_reverse:
                    reverse += 1
                    tsv_out.write(f"{read.query_name}\t{strand}\t{read.reference_start}\t{read.reference_end}\n")
                    fasta_out.write(f">{read.query_name}\n{read.query_sequence}\n")
                    reverse_bam_out.write(read)
            else:
                total += 1
                if read.is_reverse:
                    reverse += 1
                else:
                    forward += 1
                tsv_out.write(f"{read.query_name}\t{strand}\t{read.reference_start}\t{read.reference_end}\n")

    bamfile.close()
    if reverse_only:
        fasta_out.close()
        reverse_bam_out.close()
        tsv_out.close()
        sorted_bam = reverse_bam_path.replace(".bam", "_sorted.bam")
        pysam.sort("-o", sorted_bam, reverse_bam_path)
        os.replace(sorted_bam, reverse_bam_path)

    return sample_name, forward, reverse, forward + reverse if not reverse_only else reverse

def save_csv_summary(sample_stats, output_dir):
    csv_file = os.path.join(output_dir, "strand_read_summary.csv")
    with open(csv_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sample", "Forward Reads", "Reverse Reads", "Total Reads", "% Forward", "% Reverse"])
        for sample, forward, reverse, total in sample_stats:
            percent_fwd = (forward / total * 100) if total else 0
            percent_rev = (reverse / total * 100) if total else 0
            writer.writerow([sample, forward, reverse, total, f"{percent_fwd:.2f}", f"{percent_rev:.2f}"])
    print(f"\n Summary CSV saved to: {csv_file}")

def plot_summary(sample_stats, output_dir, region_start, region_end):
    samples_raw = [s[0] for s in sample_stats]
    samples = [re.sub(r'_B.*$', '', name) for name in samples_raw]
    forward = [s[1] for s in sample_stats]
    reverse = [s[2] for s in sample_stats]
    totals = [s[3] for s in sample_stats]

    fwd_percent = [f / t * 100 if t > 0 else 0 for f, t in zip(forward, totals)]
    rev_percent = [r / t * 100 if t > 0 else 0 for r, t in zip(reverse, totals)]

    x = range(len(samples))
    plt.figure(figsize=(16, 6))
    plt.bar(x, rev_percent, label='Reverse reads', color='tomato')
    plt.bar(x, fwd_percent, bottom=rev_percent, label='Forward reads', color='steelblue')

    plt.xticks(x, samples, rotation=45, ha='right')
    plt.ylabel("Percentage of Reads")
    plt.title(f"Forward and Reverse Strand Reads Spanning Region {region_start}–{region_end}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "strand_read_summary.png"))

def plot_grouped_bar(sample_stats, output_dir, region_start, region_end):
    df = pd.DataFrame(sample_stats, columns=["Sample", "Forward", "Reverse", "Total"])
    df["Sample"] = df["Sample"].apply(lambda s: re.sub(r'_B.*$', '', s))

    x = range(len(df))
    bar_width = 0.35
    plt.figure(figsize=(12, 6))
    plt.bar([i - bar_width/2 for i in x], df["Forward"], width=bar_width, label="Forward", color='steelblue')
    plt.bar([i + bar_width/2 for i in x], df["Reverse"], width=bar_width, label="Reverse", color='tomato')
    plt.xticks(x, df["Sample"], rotation=45, ha="right")
    plt.ylabel("Read Count")
    plt.title(f"Grouped Bar Chart: Reads Spanning Region {region_start}–{region_end}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "grouped_read_counts.png"))

def plot_lollipop(sample_stats, output_dir, region_start, region_end):
    df = pd.DataFrame(sample_stats, columns=["Sample", "Forward", "Reverse", "Total"])
    df["% Reverse"] = df["Reverse"] / df["Total"] * 100
    df["Sample"] = df["Sample"].apply(lambda s: re.sub(r'_B.*$', '', s))
    df = df.sort_values(by="% Reverse", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.hlines(y=df["Sample"], xmin=0, xmax=df["% Reverse"], color="gray", alpha=0.7, linewidth=2)
    plt.plot(df["% Reverse"], df["Sample"], "o", markersize=8, color="tomato")
    plt.xlabel("% Reverse Strand Reads")
    plt.title(f"Lollipop Plot: % Reverse Strand Reads (Region {region_start}–{region_end})")
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "lollipop_reverse_percent.png"))

def main():
    parser = argparse.ArgumentParser(description="Analyze strand-specific reads spanning a region.")
    parser.add_argument("bam_files", nargs="+", help="Input BAM files.")
    parser.add_argument("-o", "--output_dir", default="strand_read_output", help="Directory for output files.")
    parser.add_argument("--region", type=int, nargs=2, metavar=('START', 'END'),
                        default=(450, 633), help="Region coordinates (e.g. --region 450 633)")
    parser.add_argument("--reverse-only", action="store_true",
                        help="If set, outputs only reverse reads to TSV, FASTA, and BAM.")
    parser.add_argument("--plot-grouped", action="store_true", help="Generate grouped bar chart of read counts.")
    parser.add_argument("--plot-lollipop", action="store_true", help="Generate lollipop plot of % reverse reads.")

    args = parser.parse_args()
    region_start, region_end = args.region
    os.makedirs(args.output_dir, exist_ok=True)

    sample_stats = []
    for bam_file in args.bam_files:
        sample_name, forward, reverse, total = analyze_bam(
            bam_file, args.output_dir, region_start, region_end, reverse_only=args.reverse_only)
        sample_stats.append((sample_name, forward, reverse, total))

        print(f"\n File: {bam_file}")
        print(f"  Total spanning reads: {total}")
        if not args.reverse_only:
            print(f"    Forward: {forward} ({(forward/total)*100:.2f}%)" if total else "    Forward: 0")
        print(f"    Reverse: {reverse} ({(reverse/total)*100:.2f}%)" if total else "    Reverse: 0")

    save_csv_summary(sample_stats, args.output_dir)

    if not args.reverse_only:
        plot_summary(sample_stats, args.output_dir, region_start, region_end)
        if args.plot_grouped:
            plot_grouped_bar(sample_stats, args.output_dir, region_start, region_end)
        if args.plot_lollipop:
            plot_lollipop(sample_stats, args.output_dir, region_start, region_end)

if __name__ == "__main__":
    main()
