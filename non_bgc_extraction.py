import os
import csv
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from multiprocessing import Pool, cpu_count

# ---------- CONFIG ----------
BASE_DIR = "./"
# Easily change your input folder right here:
INPUT_FASTA_DIR = os.path.join(BASE_DIR, "whole_genome_fastas") 
metadata_file = os.path.join(BASE_DIR, "metadata_mibig.tsv")
summary_files = [
    os.path.join(BASE_DIR, "bgc_summery_mibig.csv"),
]
output_dir = os.path.join(BASE_DIR, "non_bgc_output")
os.makedirs(output_dir, exist_ok=True)

# ---------- STEP 1: Load metadata ----------
bgc_to_genome = {}  # {bgc_id: (genome_id, dataset)}
with open(metadata_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        dataset = row["dataset"].lower()
        bgc_id = row["BGC"]
        genome_id = row["genome"]
        bgc_to_genome[bgc_id] = (genome_id, dataset)

# ---------- STEP 2: Load BGC coordinates ----------
from collections import defaultdict
bgc_coords = defaultdict(list)  # {bgc_id: [(start, end), ...]}
for summary_file in summary_files:
    with open(summary_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            bgc_id = row["File Name"]
            start, end = int(row["Orig. Start"]), int(row["Orig. End"])
            if start > end:
                start, end = end, start
            bgc_coords[bgc_id].append((start, end))

# ---------- STEP 3: Group by genome ----------
genome_regions = defaultdict(list)  # {genome_id: [(start,end,bgc_id), ...]}
genome_dataset = {}
for bgc_id, (genome_id, dataset) in bgc_to_genome.items():
    if bgc_id in bgc_coords:
        for s, e in bgc_coords[bgc_id]:
            genome_regions[genome_id].append((s, e, bgc_id))
        genome_dataset[genome_id] = dataset

# ---------- STEP 4: Worker function ----------
def process_genome(genome_id):
    try:
        dataset = genome_dataset[genome_id]
        
        # Look for fasta files in the configured input folder
        genome_folder = Path(INPUT_FASTA_DIR)
        fasta_files = []
        for ext in ("*.fasta", "*.fna", "*.fa"):
            fasta_files.extend(genome_folder.glob(ext))
        
        # Try exact match first
        match_files = [f for f in fasta_files if f.stem == genome_id]
        
        # If not found, try prefix match (handles .1, .2, etc. versions)
        if not match_files:
            match_files = [f for f in fasta_files if f.stem.startswith(genome_id)]
        
        if not match_files:
            return f"[WARN] No FASTA for {genome_id}"
        
        fasta_path = match_files[0]

        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            return f"[WARN] Empty FASTA for {genome_id}"

        genome_seq = str(records[0].seq)
        regions = genome_regions[genome_id]

        # Merge overlapping intervals
        merged = []
        for start, end, _ in sorted(regions, key=lambda x: x[0]):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)

        # Build non-BGC sequence
        non_bgc_parts = []
        prev_end = 1
        for start, end in merged:
            non_bgc_parts.append(genome_seq[prev_end - 1 : start - 1])
            prev_end = end + 1
        non_bgc_parts.append(genome_seq[prev_end - 1 :])
        non_bgc_seq = "".join(non_bgc_parts)

        # Output dirs
        out_dir = Path(output_dir) / dataset
        out_dir.mkdir(parents=True, exist_ok=True)

        # Save non-BGC fasta
        nonbgc_path = out_dir / f"{genome_id}_nonBGC.fna"
        with open(nonbgc_path, "w") as out_f:
            out_f.write(f">{genome_id}_nonBGC\n")
            for i in range(0, len(non_bgc_seq), 80):
                out_f.write(non_bgc_seq[i : i + 80] + "\n")

        # Save BGC fasta
        bgc_path = out_dir / f"{genome_id}_BGCs.fna"
        with open(bgc_path, "w") as out_f:
            for idx, (start, end, bgc_id) in enumerate(sorted(regions, key=lambda x: x[0]), 1):
                seq = genome_seq[start - 1 : end]
                out_f.write(f">{genome_id}|{bgc_id}|region{idx}:{start}-{end}\n")
                for i in range(0, len(seq), 80):
                    out_f.write(seq[i : i + 80] + "\n")

        return f"[OK] {genome_id} ({dataset}) -> nonBGC:{nonbgc_path.name}, BGCs:{bgc_path.name}"

    except Exception as e:
        return f"[ERROR] {genome_id}: {e}"

# ---------- STEP 5: Run in parallel ----------
if __name__ == "__main__":
    n_cpus = max(1, cpu_count() - 2)  # use all but 2 CPUs
    with Pool(processes=n_cpus) as pool:
        results = pool.map(process_genome, genome_regions.keys())

    # Save log file
    log_file = Path(output_dir) / "processing_log.txt"
    with open(log_file, "w") as log:
        for line in results:
            print(line)
            log.write(line + "\n")

    print(f"\n[INFO] Completed. Log saved to {log_file}")