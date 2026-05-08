# operonista
Code for parsing NCBI data to visualize bacterial genomic operons :o) \
Author: Laila Phillips lphillips@utexas.edu\
References used: ChatGPT [GPT-4-turbo] / Claude [Sonnet 4.6]

An interactive command-line tool for downloading genome files from NCBI and extracting regions of interest as `.gbk` files for use with **GeneGraphics** and **clinker**.  

I use this to analyze synteny in genomic regions of bacteria, or illustrate the gene map for visualization purposes.

---

## Dependencies

| Package | Install |
|---|---|
| Biopython | `pip install biopython` |
| ncbi-datasets-cli | `conda install -c conda-forge ncbi-datasets-cli` |

---

## Usage

```bash
python genome_region_extractor.py
```

The program is fully interactive — it will prompt you for all inputs. At startup it asks whether you need to download `.gbff` files from NCBI. You can skip this if you already have them.

---

## Workflow

### Recommended workflow for multi-genome synteny (clinker)

```
1. Run program → download .gbff files into a flat directory
2. Run program → Mode C → mine that directory for your gene of interest
3. Run clinker on the output .gbk files
```

---

## Downloader

Triggered at startup with `y` when asked *"Do you need to download .gbff files from NCBI first?"*

Downloads one `.gbff` per genome accession into a flat output directory. Multi-contig genomes (draft assemblies) are handled correctly — all contig records are concatenated into a single `.gbff` file per genome.

**Accession input options:**

**Option 1 — Manual entry:** Type GCA/GCF accessions one per line.

**Option 2 — Scan a directory:** Point the tool at a folder containing files with GCA/GCF numbers in their names (e.g. from a previous NCBI download). It will automatically detect accession numbers using the standard format (`GCA_` or `GCF_` + 9 digits + optional version, e.g. `GCF_000123456.1`) and present them for confirmation before downloading.

> **Note:** NCBI GCA/GCF accession numbers always follow the format `GC[AF]_NNNNNNNNN.V` — the 9-digit portion is fixed length, so detection from filenames is reliable regardless of what else is in the filename.

---

## Extraction Modes

### Mode A — Extract by coordinates
Extracts a region from a single `.gbff` file using user-specified start/end coordinates (1-based, inclusive). If the file contains multiple contigs, you will be prompted to select one.

**Output:** `<genome>_<start>_<end>.gbk`

---

### Mode B — Extract by protein accession
Extracts the genomic region surrounding one or more protein accessions from a single `.gbff` file, with a user-defined flanking window. Version suffixes (e.g. `.1`) are stripped before matching, so either `WP_012345678` or `WP_012345678.1` will work.

**Input:** Comma-separated list of protein accessions, flank size in bp

**Output:** `<genome>_<contig>_<protein_acc>_<start>_<end>.gbk`

---

### Mode C — Search by gene/product name (multi-genome)
The main mode for synteny analysis. Searches a flat directory of `.gbff` files for a gene or product name, extracts the flanking region around every hit, and writes `.gbk` files ready for clinker.

**Expects** a flat directory where each `.gbff` is one genome — i.e. the directory produced by the built-in downloader.

**Search behavior:**
- Case-insensitive partial match
- Searches `gene`, `locus_tag`, and `product` qualifiers
- For example, querying `qseC` will match `gene="qseC"`, `product="sensor histidine kinase QseC"`, etc.

**Multiple hits per genome** are disambiguated with a numeric suffix: the first hit has none, subsequent hits are labeled `-2`, `-3`, etc.

**Input:** Path to `.gbff` directory, gene/product name, flank size in bp

**Output:** `<genome_accession>_<query>_<start>_<end>.gbk` (one file per hit)

At the end, Mode C prints a summary of genomes searched, hit counts, and the exact `clinker` command to run on your output directory.

---

## Output files

All output `.gbk` files are standard GenBank format and compatible with:
- **GeneGraphics** (https://kblin.github.io/visualize_region/) — single region visualization
- **clinker** (https://github.com/gamcil/clinker) — multi-genome synteny comparison

### Running clinker

After Mode C completes, the program prints the command for you. General form:

```bash
clinker output_dir/*.gbk -p synteny.html
```

This generates an interactive HTML figure showing gene synteny and pairwise protein identity across all extracted regions.

---

## Output filename conventions

| Mode | Filename pattern |
|---|---|
| A | `<genome>_<start>_<end>.gbk` |
| B | `<genome>_<contig>_<protein_acc>_<start>_<end>.gbk` |
| C (single hit) | `<genome_accession>_<query>_<start>_<end>.gbk` |
| C (multiple hits) | `<genome_accession>_<query>-2_<start>_<end>.gbk` |

---

## Notes

- `.gbff` files from NCBI are often multi-record (one record per contig). All modes handle this correctly — the original script only read the first contig, which caused silent data loss for draft assemblies.
- Mode B strips version suffixes before matching protein accessions, fixing a bug where e.g. `WP_012345678.1` would not match if entered without the version number or vice versa.
- Mode C is non-recursive by design — it only reads `.gbff` files directly inside the specified directory, not subdirectories. This avoids accidental parsing of unrelated files in nested genome directory trees.
