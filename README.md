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

# What next?
The file output(s) containing protein sequences and easily-identifiable fasta headers should be in an acceptable format for most tree-building pipelines. Here is what I do with mine, for example:
### 1. Align protein sequences
Use an alignment tool to align all of the protein sequences pulled from operonista. I usually use command-line MUSCLE in a conda environment, but there is also a UI version:\
https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein\ \
\
**CAUTION:** operonista is pulling sequences based on the annotation provided by Operon-mapper, and I cannot guarantee these annotations are fully accurate and/or that the protein-coding genes are homologous. This may leave you with a chaotic phylogeny if there is little similarity between the genes you pull. In these situations, you may consider building that chaotic tree first, and using that to hone in on particular branches of interest (e.g., build another tree with that subset).
### 2. Trim alignments (optional)
Trimming alignments can be useful to remove gaps in your alignments that don't offer much resolution for your phylogeny. I usually use ClipKit:\
https://github.com/JLSteenwyk/ClipKIT
### 3. Build tree
Once you have your fasta file containing aligned and trimmed protein sequences, you're ready for the big guns. I usually use IQTree to build my phylogenies, with options like ModelFinder and bootstrapping:
http://www.iqtree.org/
