import re

# File paths
genome_ids_file = "/Users/Snacks/Grad_School/Moran_Lab/Projects/Orbaceae/whole_genome_sequencing/operon-mapper/PapC/genome_IDs.txt"
fasta_file = "/Users/Snacks/Grad_School/Moran_Lab/Projects/Orbaceae/whole_genome_sequencing/operon-mapper/PapC/combined_PapC_sequences.fasta"
output_fasta_file = "/Users/Snacks/Grad_School/Moran_Lab/Projects/Orbaceae/whole_genome_sequencing/operon-mapper/PapC/combined_PapC_sequences_renamed.fasta"

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

