import re

# File paths
genome_ids_file = "genome_IDs.txt" # change to filepath for your tab-separated list of genome IDs and genome names
fasta_file = "combined_pilus_sequences.fasta" # change to filepath for the target fasta file you want to rename
output_fasta_file = "combined_pilus_sequences_renamed.fasta" # specify name for your output file

# Read genome_IDs.txt into a dictionary
prefix_to_genome = {}
with open(genome_ids_file, "r") as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 3:
            _, prefix, genome = parts
            prefix_to_genome[prefix] = genome

# Process the FASTA file and update headers
with open(fasta_file, "r") as infile, open(output_fasta_file, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            match = re.match(r">(\S+?)_(\S+)", line)  # Extract prefix before "_"
            if match:
                protein_id, rest = match.groups()
                if protein_id in prefix_to_genome:
                    new_header = f">{prefix_to_genome[protein_id]}_{rest}\n"
                    outfile.write(new_header)
                else:
                    outfile.write(line)  # Keep unchanged if no match found
            else:
                outfile.write(line)  # Keep unchanged if format is unexpected
        else:
            outfile.write(line)  # Write sequence lines unchanged

print(f"Updated FASTA file saved as: {output_fasta_file}")

