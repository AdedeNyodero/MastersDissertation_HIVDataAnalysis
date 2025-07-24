import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os

def parse_blast_file(filenames):
    blast_data = []
    for file in filenames:
        with open(file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                try:
                    read_id = parts[0]
                    reference = parts[1]
                    ref_start = int(parts[2])
                    ref_end = int(parts[3])
                    query_start = int(parts[4])
                    query_end = int(parts[5])
                    ref_seq = parts[6]
                    query_seq = parts[7]
                except ValueError:
                    continue
                blast_data.append({
                    'read_id': read_id,
                    'reference': reference,
                    'ref_start': ref_start,
                    'ref_end': ref_end,
                    'query_start': query_start,
                    'query_end': query_end,
                    'ref_seq': ref_seq,
                    'query_seq': query_seq
                })
    return blast_data

def parse_lengths_file(length_file):
    lengths = {}
    with open(length_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 2:
                continue
            read_id, length_str = parts[0], parts[1]
            try:
                lengths[read_id] = int(length_str)
            except ValueError:
                continue
    return lengths

def infer_strand(ref_start, ref_end):
    return '+' if ref_end > ref_start else '-'

def infer_overhang(q_start, q_end, read_length, strand):
    if strand == '+':
        if q_start > 1:
            return "5' overhang"
        elif read_length > q_end:
            return "3' overhang"
    else:
        if q_end < read_length:
            return "5' overhang"
        elif q_start > 1:
            return "3' overhang"
    return "no overhang"

def bin_position(pos, bin_size):
    return (pos // bin_size) * bin_size

def count_reads(blast_data, lengths, bin_size=100):
    counts = defaultdict(lambda: defaultdict(int))
    for row in blast_data:
        read_id = row['read_id']
        if read_id not in lengths:
            continue
        read_length = lengths[read_id]
        strand = infer_strand(row['ref_start'], row['ref_end'])
        overhang = infer_overhang(row['query_start'], row['query_end'], read_length, strand)
        bin_pos = bin_position(row['ref_start'], bin_size)
        counts[(strand, overhang)][bin_pos] += 1
    return counts

def plot_counts(counts, bin_size, output_file):
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), sharey=True)

    color_map = {"5' overhang": "#1C5AE9", "3' overhang": "#E91212"}

    labels = [
        ('+', "5' overhang", axs[0, 0], "Forward templates with 5' overhang"),
        ('+', "3' overhang", axs[0, 1], "Forward templates with 3' overhang"),
        ('-', "5' overhang", axs[1, 0], "Reverse templates with 5' overhang"),
        ('-', "3' overhang", axs[1, 1], "Reverse templates with 3' overhang")
    ]

    for strand, overhang, ax, title in labels:
        bin_counts = counts.get((strand, overhang), {})
        bins = sorted(bin_counts)
        values = [bin_counts[b] for b in bins]
        ax.bar(bins, values, width=bin_size * 0.8, align='edge', color=color_map.get(overhang, 'gray'))
        ax.set_title(title)
        ax.set_xlabel("Genomic position (binned)")
        #ax.grid(True, linestyle='--', alpha=0.5)
        if bins:
            ax.set_xlim(min(bins), max(bins) + bin_size)

    axs[0, 0].set_ylabel("Read counts")
    axs[1, 0].set_ylabel("Read counts")
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Full overhang plot saved to {output_file}")

def plot_strand_overhang(counts, bin_size, strand, overhang, color, output_file):
    bin_counts = counts.get((strand, overhang), {})
    bins = sorted(bin_counts)
    values = [bin_counts[b] for b in bins]

    if not bins:
        print(f"No reads found for {strand} strand with {overhang}")
        return

    plt.figure(figsize=(12, 6))
    plt.bar(bins, values, width=bin_size * 0.8, color=color)
    strand_label = "Forward" if strand == '+' else "Reverse"
    plt.title(f"{strand_label} strand reads with {overhang}")
    plt.xlabel("Genomic position (binned)")
    plt.ylabel("Read counts")
    plt.xticks(ticks=bins, labels=[f"{b}-{b + bin_size - 1}" for b in bins], rotation=90, ha='right')
    #plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"{strand_label} strand {overhang} plot saved to {output_file}")

def plot_reverse_strand_only(counts, bin_size, output_file):
    overhangs = ["5' overhang", "3' overhang"]
    colors = {"5' overhang": "#1C5AE9", "3' overhang": "#E91212"}
    bin_set = set()

    for overhang in overhangs:
        bin_set.update(counts.get(('-', overhang), {}).keys())
    bins = sorted(bin_set)

    if not bins:
        print("No reverse strand reads to plot.")
        return

    data = {o: [counts.get(('-', o), {}).get(b, 0) for b in bins] for o in overhangs}
    x = list(range(len(bins)))
    bar_width = 0.4

    plt.figure(figsize=(12, 6))
    for i, overhang in enumerate(overhangs):
        shift = -bar_width/2 if i == 0 else bar_width/2
        plt.bar(
            [xi + shift for xi in x],
            data[overhang],
            width=bar_width,
            label=overhang,
            color=colors[overhang]
        )

    plt.xticks(ticks=x, labels=[f"{b}-{b + bin_size - 1}" for b in bins], rotation=90, ha='right')
    plt.xlabel("Genomic position (binned)")
    plt.ylabel("Read counts")
    plt.title("Reverse Strand Read Distribution by Overhang Type")
    plt.legend()
    plt.tight_layout()

    outname, ext = os.path.splitext(output_file)
    reverse_output = f"{outname}_reverse{ext}"
    plt.savefig(reverse_output)
    print(f"Reverse strand plot saved to {reverse_output}")

def plot_forward_strand_separate(counts, bin_size, output_file_base):
    color_map = {"5' overhang": "#1C5AE9", "3' overhang": "#E91212"}
    for overhang in ["5' overhang", "3' overhang"]:
        outname, ext = os.path.splitext(output_file_base)
        suffix = '5prime' if overhang == "5' overhang" else '3prime'
        output = f"{outname}_forward_{suffix}{ext}"
        plot_strand_overhang(counts, bin_size, '+', overhang, color_map[overhang], output)

def main():
    parser = argparse.ArgumentParser(description="Plot HIV read overhangs with separate strand visualizations.")
    parser.add_argument('--blast', nargs='+', required=True, help='BLAST tabular output files (no headers).')
    parser.add_argument('--lengths', required=True, help='Read lengths file (CSV, no header).')
    parser.add_argument('--output', required=True, help='Base name for output plots (e.g. overhangs.png)')
    parser.add_argument('--binsize', type=int, default=100, help='Bin size for genomic position')

    args = parser.parse_args()

    blast_data = parse_blast_file(args.blast)
    lengths = parse_lengths_file(args.lengths)
    counts = count_reads(blast_data, lengths, args.binsize)

    plot_counts(counts, args.binsize, args.output)
    plot_reverse_strand_only(counts, args.binsize, args.output)
    plot_forward_strand_separate(counts, args.binsize, args.output)

if __name__ == '__main__':
    main()
