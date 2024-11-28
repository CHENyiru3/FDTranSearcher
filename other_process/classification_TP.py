import os
from Bio import SeqIO
from collections import defaultdict

def organize_fasta_by_classification_combined(fasta_file, classification_file, output_base_dir):
    """
    Organizes protein sequences from a FASTA file based on classifications,
    combines sequences of the same subclass into a single FASTA file,
    and stores them in a structured directory hierarchy.

    Args:
        fasta_file (str): Path to the FASTA file containing protein sequences.
        classification_file (str): Path to the file containing protein classifications.
        output_base_dir (str): Base directory to create the output directory structure.
    """

    # 1. Read Classifications
    classifications = {}
    with open(classification_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            rexdb_id = parts[0]
            classification = parts[1:]
            classifications[rexdb_id] = classification

    # 2. Parse FASTA and Organize Sequences (Combined)
    organized_sequences = defaultdict(list)  # Use a dictionary to group sequences by classification

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        try:
            rexdb_id = header.split('__')[-1]
        except IndexError:
            print(f"Warning: Could not extract REXdb_ID from header: {header}. Skipping.")
            continue

        if rexdb_id in classifications:
            classification_levels = classifications[rexdb_id]
            # Use the full classification path as the key to group sequences
            classification_key = os.path.join(*classification_levels)
            organized_sequences[classification_key].append(record)
        else:
            print(f"Warning: No classification found for {rexdb_id}. Skipping.")

    # 3. Write Combined FASTA Files
    for classification_key, sequences in organized_sequences.items():
        output_dir = os.path.join(output_base_dir, *classification_key.split(os.sep)[:-1]) # Get directory path without filename
        os.makedirs(output_dir, exist_ok=True)
        output_filename = os.path.join(output_dir, f"{classification_key.split(os.sep)[-1]}.fasta")  # Use last level of classification as filename

        with open(output_filename, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

# Example Usage (same as before)
fasta_file = "/mnt/volume1/2023SRTP/data/maize/Viridiplantae_v4.0.fasta"
classification_file = "/mnt/volume1/2023SRTP/data/maize/Viridiplantae_v4.classification"
output_base_dir = "/mnt/volume1/2023SRTP/data/maize/protein_family"

organize_fasta_by_classification_combined(fasta_file, classification_file, output_base_dir)
print("Protein sequences organized and combined successfully!")