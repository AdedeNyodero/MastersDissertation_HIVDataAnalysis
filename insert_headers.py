import pandas as pd
import os
from glob import glob

# Folder where the .tsv files are located
input_folder = "." 

# Defining headers
headers = [
    "Read name",
    "Reference name",
    "Reference mapping start",
    "Reference mapping end",
    "Query start",
    "Query End",
    "Reference alignment",
    "Query alignment"
]

# Get all .tsv files
tsv_files = glob(os.path.join(input_folder, "*.tsv"))

for file in tsv_files:
    if "with_headers" in file:
        continue  # Skip already processed files

    try:
        # Read the file with no header
        df = pd.read_csv(file, sep="\t", header=None)

        # Add the headers
        df.columns = headers

        # Save to new TSV with header
        output_file = file.replace(".tsv", "_with_headers.tsv")
        df.to_csv(output_file, sep="\t", index=False)

        print(f" Processed: {os.path.basename(output_file)}")

    except Exception as e:
        print(f" Failed to process {file}: {e}")
