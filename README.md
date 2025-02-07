# operonista
Code for parsing operon data from bacterial genomes
Author: Laila Phillips lphillips@utexas.edu
References used: ChatGPT [GPT-4-turbo]

# Dependencies
Python

# Input
Input files are the output from Operon-mapper (https://biocomputo.ibt.unam.mx/operon_mapper/)

# Directory structure:
Note: The only filename that needs to be explicitly specified in the script is the root directory.

root_directory/  # The directory specified in the script
├── abc123/
│   ├── list_of_operons_abc123.txt
│   ├── predicted_protein_sequences_abc123.txt
├── another_genome/
│   ├── list_of_operons_genome2.txt
│   ├── predicted_protein_sequences_genome2.txt
├── random_folder/
│   ├── list_of_operons_xyz789.txt
│   ├── predicted_protein_sequences_xyz789.txt
