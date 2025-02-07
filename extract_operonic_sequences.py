import csv
import os
import glob


def extract_sequences_from_directory(root_dir, output_fasta, output_verification, keyword="pilus"):
    """
    Searches through a directory structure for operon and sequence files, extracts sequences
    based on matching protein IDs, and compiles results into a single FASTA file.

    Parameters:
    - root_dir: Directory to search through (contains subdirectories with genome files).
    - output_fasta: Path to save the combined FASTA file.
    - output_verification: Path to save the verification file with ProteinID and Function.
    - keyword: Keyword to search for in the Function column (default: "pilus").
    """

    all_sequences = []  # Store extracted FASTA entries
    all_annotations = []  # Store ProteinID and Function for verification

    # Recursively find all matching operon and sequence files in subdirectories
    for subdir, _, _ in os.walk(root_dir):
        # Debugging: Print which directory is being checked
        print(f"Checking directory: {subdir}")

        operon_files = glob.glob(os.path.join(subdir, "list_of_operons_*"))
        sequence_files = glob.glob(os.path.join(subdir, "predicted_protein_sequences_*"))

        # Debugging: Print found files
        print(f"Found operon files: {operon_files}")
        print(f"Found sequence files: {sequence_files}")

        # Ensure both required files exist in the subdirectory
        if not operon_files or not sequence_files:
            print(f"Skipping {subdir} (missing required files)")
            continue  # Skip if files are missing

        # Assume one operon file and one sequence file per genome
        operon_file = operon_files[0]
        sequence_file = sequence_files[0]

        print(f"Processing: {operon_file} and {sequence_file}")

        # Step 1: Extract matching Protein IDs and their annotations
        matching_protein_ids = {}
        with open(operon_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # Skip header

            for row in reader:
                #print(f"Row: {row}")  # Debugging
                if len(row) >= 8:
                    protein_id, function = row[1], row[7]
                    if keyword.lower() in function.lower():
                        matching_protein_ids[protein_id] = function  # Store function as well

        # Step 2: Extract sequences for matching IDs
        with open(sequence_file, "r") as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if lines[i].startswith(">"):
                    protein_id = lines[i][1:].strip()
                    if protein_id in matching_protein_ids:
                        all_sequences.append(lines[i])  # Header
                        all_sequences.append(lines[i + 1])  # Sequence
                        all_annotations.append(f"{protein_id}\t{matching_protein_ids[protein_id]}\n")

    # Step 3: Write results to output files
    with open(output_fasta, "w") as f:
        f.writelines(all_sequences)

    with open(output_verification, "w") as f:
        f.write("ProteinID\tFunction\n")  # Add header
        f.writelines(all_annotations)

    print(f"Extracted {len(all_sequences) // 2} sequences containing '{keyword}'.")
    print(f"FASTA saved to: {output_fasta}")
    print(f"Verification file saved to: {output_verification}")


# Example usage:
if __name__ == "__main__":
    root_directory = "genome_data"  # Change to your main directory containing subdirectories
    output_fasta_file = "combined_pilus_sequences.fasta"
    output_verification_file = "verification_pilus_sequences.txt"
    keyword = "pilus"  # Modify if needed

    extract_sequences_from_directory(root_directory, output_fasta_file, output_verification_file, keyword)