# operonista
Code for parsing Operon-mapper output from bacterial genomes :o) \
Author: Laila Phillips lphillips@utexas.edu\
References used: ChatGPT [GPT-4-turbo]

## Dependencies
1. Operon-mapper [https://biocomputo.ibt.unam.mx/operon_mapper/]
2. Python

## Purpose
This code is intended to pull out protein sequences from the list_of_operons_* and predicted_protein_sequences_* output files matching keywords you specify in the extract_operonic_sequences.py file. This is useful for phylogenetic analysis of certain genes.

## Inputs
1. Input files for grabbing sequences matching keyword are the outputs from Operon-mapper
2. Input file for renaming fasta headers is a tab-separated file containing the following:
  >"operon-mapperID" (a number generated for each Operon-mapper job)\
  >"proteinIDprefix" (an alphabetic code specific to each genome that is a prefix for gene ids within Operon-mapper output files)\
  >"genome" (the name you want to swap in for your fasta headers)
### Directory structure:
Note: The only filename that needs to be explicitly specified in the script is the root directory.

root_directory/  # The directory specified in the script\
├── abc123/\
│   ├── list_of_operons_abc123.txt\
│   ├── predicted_protein_sequences_abc123.txt\
├── another_genome/\
│   ├── list_of_operons_genome2.txt\
│   ├── predicted_protein_sequences_genome2.txt\
├── random_folder/\
│   ├── list_of_operons_xyz789.txt\
│   ├── predicted_protein_sequences_xyz789.txt\
