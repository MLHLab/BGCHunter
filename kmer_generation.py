import os
import csv
import itertools
import re  # Added for splitting sequences
from multiprocessing import Pool, cpu_count
import numpy as np

# Input and output
folder_path = "./mibig_fasta_outputs"
output_file = "mibig_bgc_6mer.csv"
K = 6

# Generate all possible kmers (AAAAAA...TTTTTT)
bases = ["A", "C", "G", "T"]
all_kmers = ["".join(p) for p in itertools.product(bases, repeat=K)]
kmer_index = {kmer: i for i, kmer in enumerate(all_kmers)}

def count_kmers_in_file(file_path):
    """Count overlapping 6-mers in a single .fna file."""
    sequence_lines = []
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                # Just append, do not replace "N" here
                sequence_lines.append(line.strip().upper()) 
    
    full_sequence = "".join(sequence_lines)
    counts = np.zeros(len(all_kmers), dtype=int)

    # Split the sequence by 'N' (or any non-ACGT character) to avoid phantom kmers
    valid_chunks = re.split(r'[^ACGT]', full_sequence)

    # Run your original counting logic, but on each valid chunk separately
    for chunk in valid_chunks:
        if len(chunk) < K:
            continue  # Skip chunks that are too small to contain a k-mer
            
        for i in range(len(chunk) - K + 1):
            kmer = chunk[i:i+K]
            if kmer in kmer_index:  # only A/C/G/T kmers
                counts[kmer_index[kmer]] += 1
    
    return os.path.basename(file_path), counts

if __name__ == "__main__":
    fna_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".fna")]

    # Use all available CPUs (or limit if needed)
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(count_kmers_in_file, fna_files)

    # Write CSV
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File"] + all_kmers)
        for filename, counts in results:
            writer.writerow([filename] + counts.tolist())

    print(f"✅ Parallel {K}-mer counts saved to {output_file}")