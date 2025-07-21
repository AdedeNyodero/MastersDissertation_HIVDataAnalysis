import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import gzip

def parse_fastq(fastq_file):
    """Reads sequences from gzipped or plain FASTQ file into a dictionary."""
    read_seqs = {}
    open_func = gzip.open if fastq_file.endswith(".gz") else open
    with open_func(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_seqs[record.id] = str(record.seq)
    return read_seqs

def determine_strand(flag):
    return 'reverse' if (flag & 16) != 0 else 'forward'

def analyze_overhangs(bam_file, fastq_file, bin_size=100, genome_length=9719):
    bam = pysam.AlignmentFile(bam_file, "rb")
    fastq_reads = parse_fastq(fastq_file)

    overhang_data = []

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        read_name = read.query_name
        strand = determine_strand(read.flag)
        ref_start = read.reference_start
        ref_end = read.reference_end
        aligned_seq = read.query_sequence
        full_seq = fastq_reads.get(read_name)

        if not full_seq or not aligned_seq:
            continue

        overhang_5p = 0
        overhang_3p = 0
        overhang_type = "none"

        # Compare aligned sequence to full sequence
        if strand == 'forward':
            start_index = full_seq.find(aligned_seq)
            if start_index != -1:
                overhang_5p = start_index
                overhang_3p = len(full_seq) - (start_index + len(aligned_seq))
        else:  # reverse strand
            rev_full = str(Seq(full_seq).reverse_complement())
            start_index = rev_full.find(aligned_seq)
            if start_index != -1:
                overhang_5p = start_index
                overhang_3p = len(rev_full) - (start_index + len(aligned_seq))

        if overhang_5p > 0 and overhang_3p > 0:
            overhang_type = "both"
        elif overhang_5p > 0:
            overhang_type = "5p"
        elif overhang_3p > 0:
            overhang_type = "3p"

        overhang_data.append({
            'read_name': read_name,
            'strand': strand,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'overhang_5p': overhang_5p,
            'overhang_3p': overhang_3p,
            'overhang_type': overhang_type,
        })

    bam.close()
    return pd.DataFrame(overhang_data)

def plot_overhangs(df, bin_size, genome_length, output_prefix):
    df['bin'] = (df['ref_start'] // bin_size) * bin_size
    palette = {'5p': 'blue', '3p': 'orange', 'both': 'purple'}

    for strand in ['forward', 'reverse']:
        strand_df = df[df['strand'] == strand]
        plot_df = strand_df.groupby(['bin', 'overhang_type']).size().reset_index(name='count')

        plt.figure(figsize=(14, 5))
        sns.barplot(data=plot_df, x='bin', y='count', hue='overhang_type', palette=palette)
        plt.title(f"Read Overhangs on {strand.capitalize()} Strand")
        plt.xlabel('HIV Genome Position (binned)')
        plt.ylabel('Number of Reads with Overhangs')
        plt.xticks(rotation=45)
        plt.legend(title='Overhang Type')
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_{strand}_overhangs.png", dpi=300)
        plt.close()

def plot_overhang_distributions(df, bin_size, genome_length, output_prefix):
    df['bin'] = (df['ref_start'] // bin_size) * bin_size
    df['bin_label'] = df['bin'].astype(str) + '-' + (df['bin'] + bin_size - 1).astype(str)

    fig, axs = plt.subplots(2, 2, figsize=(18, 10), sharex=True, sharey=True)

    for i, (strand, overhang_end, color) in enumerate([
        ('forward', '5p', 'blue'),
        ('forward', '3p', 'red'),
        ('reverse', '5p', 'blue'),
        ('reverse', '3p', 'red')
    ]):
        ax = axs[i // 2, i % 2]

        # Filter for strand
        strand_df = df[df['strand'] == strand]
        total_per_bin = strand_df.groupby('bin_label').size()

        # Filter reads with overhang on specific end
        overhang_col = 'overhang_5p' if overhang_end == '5p' else 'overhang_3p'
        filtered_df = strand_df[strand_df[overhang_col] > 0]
        overhang_per_bin = filtered_df.groupby('bin_label').size()

        # Proportion
        proportion = (overhang_per_bin / total_per_bin).fillna(0).sort_index()

        ax.plot(proportion.index, proportion.values, color=color, marker='x', linestyle='-')
        ax.set_ylim(0, 1)
        ax.set_xticks(range(len(proportion.index)))
        ax.set_xticklabels(proportion.index, rotation=90, fontsize=7)
        ax.set_title(f"Distribution of {strand} read % with {overhang_end} overhang")
        ax.set_ylabel("Proportion of reads")

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_overhang_distributions.png", dpi=300)
    plt.close()
    print(f"\nSaved overhang distribution plot: {output_prefix}_overhang_distributions.png")

def write_summary(df, output_prefix):
    summary = df.groupby(['strand', 'overhang_type']).size().reset_index(name='read_count')
    summary.to_csv(f"{output_prefix}_summary.csv", index=False)
    print("\nSummary statistics:")
    print(summary)

def write_detailed_csv(df, output_prefix):
    df.to_csv(f"{output_prefix}_overhangs.csv", index=False)
    print(f"\nDetailed per-read data written to: {output_prefix}_overhangs.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect HIV read overhangs and visualize")
    parser.add_argument("--fastq", required=True, help="Original FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--bam", required=True, help="Filtered BAM file (HIV-only reads)")
    parser.add_argument("--genome_length", type=int, default=9719, help="Length of HIV genome")
    parser.add_argument("--bin_size", type=int, default=100, help="Bin size for x-axis")
    parser.add_argument("--output_prefix", default="hiv_overhang", help="Prefix for output files")

    args = parser.parse_args()

    df = analyze_overhangs(args.bam, args.fastq, args.bin_size, args.genome_length)
    write_detailed_csv(df, args.output_prefix)
    write_summary(df, args.output_prefix)
    plot_overhangs(df, args.bin_size, args.genome_length, args.output_prefix)
    plot_overhang_distributions(df, args.bin_size, args.genome_length, args.output_prefix)
