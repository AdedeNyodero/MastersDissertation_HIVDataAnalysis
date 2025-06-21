#!/usr/bin/env python3

import os
import re
import matplotlib.pyplot as plt
import numpy as np


def parse_summary_file(path):
    with open(path) as f:
        lines = f.readlines()
    sample = os.path.basename(os.path.dirname(path))
    forward_line = next(line for line in lines if "Forward strand" in line)
    reverse_line = next(line for line in lines if "Reverse strand" in line)

    forward_pct = float(re.search(r"\((.*?)%\)", forward_line).group(1))
    reverse_pct = float(re.search(r"\((.*?)%\)", reverse_line).group(1))

    return sample, forward_pct, reverse_pct


def collect_data(summary_dir):
    samples = []
    forward_pcts = []
    reverse_pcts = []

    for root, dirs, files in os.walk(summary_dir):
        for file in files:
            if file == "strand_summary.txt":
                sample, fwd, rev = parse_summary_file(os.path.join(root, file))
                samples.append(sample)
                forward_pcts.append(fwd)
                reverse_pcts.append(rev)

    return samples, forward_pcts, reverse_pcts


def plot_strand_proportions(samples, forward_pcts, reverse_pcts, output_file):
    x = np.arange(len(samples))
    width = 0.5

    fig, ax = plt.subplots(figsize=(16, 6))
    ax.bar(x, reverse_pcts, width, label='Reverse reads', color='tomato')
    ax.bar(x, forward_pcts, width, bottom=reverse_pcts, label='Forward reads', color='steelblue')

    for i in range(len(samples)):
        ax.text(i, reverse_pcts[i]/2, f"{reverse_pcts[i]:.1f}%", ha='center', va='center', color='white', fontsize=10, fontweight='bold')
        ax.text(i, reverse_pcts[i] + forward_pcts[i]/2, f"{forward_pcts[i]:.1f}%", ha='center', va='center', color='white', fontsize=10, fontweight='bold')

    ax.set_ylabel('Percentage of Reads')
    ax.set_title('Proportion of Reverse vs Forward Reads per Dataset')
    ax.set_xlabel('Datasets')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.set_ylim(0, 110)
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Plot reverse vs forward strand proportions from summary files.")
    parser.add_argument("--summary_root", required=True, help="Top-level directory containing sample output folders")
    parser.add_argument("--output", default="strand_proportions.png", help="Output plot filename")
    args = parser.parse_args()

    samples, forward_pcts, reverse_pcts = collect_data(args.summary_root)
    plot_strand_proportions(samples, forward_pcts, reverse_pcts, args.output)
    print(f" Saved plot to {args.output}")


if __name__ == "__main__":
    main()
