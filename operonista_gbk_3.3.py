#!/usr/bin/env python3
"""
=======================
OPERONISTA
VERSION: 3.3
LAST UPDATED: 5/17/26
=======================
Integrated tool for:
  - Downloading .gbff files from NCBI by GCA/GCF accession
  - Extracting genomic regions by coordinate, protein accession, or gene/product name
  - Output .gbk files compatible with GeneGraphics and clinker

Requires:
  - biopython         (pip install biopython)
  - ncbi-datasets-cli (conda install -c conda-forge ncbi-datasets-cli)
"""

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import copy
import os
import sys
import subprocess
import shutil
import zipfile
import glob
import re


# ─────────────────────────────────────────────
# NCBI DOWNLOADER
# ─────────────────────────────────────────────

def check_datasets_cli():
    """Confirm ncbi-datasets-cli is available."""
    if shutil.which("datasets") is None:
        print(
            "\n[ERROR] 'datasets' command not found.\n"
            "Install it with: conda install -c conda-forge ncbi-datasets-cli\n"
        )
        sys.exit(1)

def download_gbff(accessions, download_dir):
    """
    Download .gbff files for a list of GCA/GCF accessions using ncbi datasets CLI.
    Each genome's .gbff is saved as <accession>.gbff in download_dir.
    """
    check_datasets_cli()
    os.makedirs(download_dir, exist_ok=True)

    print(f"\nDownloading {len(accessions)} genome(s) to: {download_dir}")

    for acc in accessions:
        acc = acc.strip()
        if not acc:
            continue

        zip_path = os.path.join(download_dir, f"{acc}.zip")
        print(f"  Fetching {acc} ...", end=" ", flush=True)

        result = subprocess.run(
            [
                "datasets", "download", "genome", "accession", acc,
                "--include", "gbff",
                "--filename", zip_path,
                "--no-progressbar",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print(f"FAILED\n    {result.stderr.strip()}")
            continue

        # Unzip and relocate the .gbff file
        try:
            with zipfile.ZipFile(zip_path, "r") as zf:
                gbff_files = [f for f in zf.namelist() if f.endswith(".gbff")]
                if not gbff_files:
                    print("FAILED (no .gbff found in archive)")
                    continue
                # There may be multiple .gbff files (one per contig batch); concatenate them
                out_path = os.path.join(download_dir, f"{acc}.gbff")
                with open(out_path, "wb") as out_f:
                    for gbff_name in gbff_files:
                        with zf.open(gbff_name) as src:
                            shutil.copyfileobj(src, out_f)
            os.remove(zip_path)
            print(f"OK → {acc}.gbff")
        except zipfile.BadZipFile:
            print("FAILED (bad zip file)")

    print("\nDownload complete.\n")


# GCA/GCF accession pattern: GCA_ or GCF_ + exactly 9 digits + optional .version
_ACCESSION_RE = re.compile(r"(GC[AF]_\d{9}(?:\.\d+)?)", re.IGNORECASE)


def scan_filenames_for_accessions(directory):
    """
    Scan all filenames (not contents) in a directory for GCA/GCF accession numbers.
    Returns a sorted, deduplicated list of accessions found.
    Non-recursive — only looks at the immediate directory listing.
    """
    found = {}  # accession_without_version -> full accession (keep longest/last seen)
    try:
        entries = os.listdir(directory)
    except FileNotFoundError:
        print(f"  [ERROR] Directory not found: {directory}")
        return []

    for name in entries:
        for match in _ACCESSION_RE.finditer(name):
            acc = match.group(1)
            bare = acc.split(".")[0].upper()
            # If we see the same base accession twice with different versions, keep
            # the one with a version number (more specific)
            existing = found.get(bare)
            if existing is None or "." not in existing:
                found[bare] = acc

    return sorted(found.values())


def run_downloader():
    """Interactive downloader entry point."""
    print("\n─── NCBI .gbff Downloader ───")
    print("\nHow would you like to provide accession numbers?")
    print("  1 — Type them in manually")
    print("  2 — Scan filenames in a directory (auto-detect GCA/GCF numbers)")
    source = input("Choose [1/2]: ").strip()

    accessions = []

    if source == "2":
        scan_dir = input("  Directory to scan: ").strip()
        accessions = scan_filenames_for_accessions(scan_dir)
        if not accessions:
            print("  No GCA/GCF accessions found in filenames. Exiting downloader.")
            return
        print(f"\n  Found {len(accessions)} accession(s):")
        for a in accessions:
            print(f"    {a}")
        confirm = input("\n  Proceed with these? [Y/n]: ").strip().lower()
        if confirm == "n":
            print("  Cancelled.")
            return

    else:
        print("Enter GCA/GCF accessions one per line.")
        print("When done, enter a blank line or type 'done'.\n")
        while True:
            line = input("  Accession: ").strip()
            if line.lower() in ("", "done"):
                break
            if _ACCESSION_RE.search(line):
                accessions.append(line.strip())
            else:
                print(f"  [WARN] '{line}' doesn't look like a GCA/GCF accession — skipping.")

    if not accessions:
        print("No accessions to download. Skipping.")
        return

    download_dir = input(
        "\nDirectory to save .gbff files (will be created if needed): "
    ).strip()
    download_gbff(accessions, download_dir)


# ─────────────────────────────────────────────
# GENBANK PARSING UTILITIES
# ─────────────────────────────────────────────

def load_all_records(gbff_path):
    """
    Parse all records from a .gbff file (multi-record / multi-contig safe).
    Returns list of SeqRecord objects.
    """
    records = list(SeqIO.parse(gbff_path, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in: {gbff_path}")
    return records


def features_overlap(f_start, f_end, region_start, region_end):
    return (f_end > region_start) and (f_start < region_end)


def shift_location(loc, shift, region_len):
    """Shift and clip FeatureLocation or CompoundLocation to region bounds."""
    if isinstance(loc, CompoundLocation):
        parts = [shift_location(p, shift, region_len) for p in loc.parts]
        parts = [p for p in parts if p is not None]
        if not parts:
            return None
        if len(parts) == 1:
            # CompoundLocation requires ≥2 parts; downgrade to SimpleLocation
            return parts[0]
        return CompoundLocation(parts, loc.operator)

    start = int(loc.start) + shift
    end = int(loc.end) + shift

    if end <= 0 or start >= region_len:
        return None
    start = max(0, start)
    end = min(region_len, end)

    return FeatureLocation(start, end, strand=loc.strand)


def extract_region(record, start, end):
    """
    Extract a region from a SeqRecord (1-based, inclusive coords).
    Returns a new SeqRecord with correctly shifted features.
    """
    region_start = start - 1  # convert to 0-based
    region_end = end
    region_len = region_end - region_start

    sub = record[region_start:region_end]
    new_features = []

    for feat in record.features:
        f_start = int(feat.location.start)
        f_end = int(feat.location.end)

        if not features_overlap(f_start, f_end, region_start, region_end):
            continue

        new_feat = copy.deepcopy(feat)
        new_loc = shift_location(feat.location, -region_start, region_len)
        if new_loc is None:
            continue

        new_feat.location = new_loc
        new_features.append(new_feat)

    sub.features = new_features
    return sub


def safe_record_id(record):
    """Return a filesystem-safe record ID (strip version suffix for filenames)."""
    return record.id.replace(".", "_")


# ─────────────────────────────────────────────
# MODE A — COORDINATE EXTRACTION
# ─────────────────────────────────────────────

def run_mode_a(records, accession_label, output_dir):
    print("\n─── Mode A: Extract by Coordinates ───")

    if len(records) > 1:
        print(f"  Note: {len(records)} contigs in this file.")
        for i, r in enumerate(records):
            print(f"    [{i}] {r.id}  ({len(r)} bp)")
        idx = int(input("  Select contig index: ").strip())
        record = records[idx]
    else:
        record = records[0]

    start = int(input("  Start coordinate: ").strip())
    end   = int(input("  End coordinate:   ").strip())
    if start > end:
        start, end = end, start

    sub = extract_region(record, start, end)
    outname = os.path.join(output_dir, f"{accession_label}_{start}_{end}.gbk")
    SeqIO.write(sub, outname, "genbank")
    print(f"  ✔ Written: {outname}")


# ─────────────────────────────────────────────
# MODE B — PROTEIN ACCESSION EXTRACTION (fixed)
# ─────────────────────────────────────────────

def run_mode_b(records, accession_label, output_dir):
    print("\n─── Mode B: Extract by Protein Accession ───")

    prot_input = input(
        "  Enter protein accession(s), comma-separated: "
    ).strip()
    prot_list = [p.strip() for p in prot_input.split(",")]
    window = int(input("  Flank size (bp): ").strip())

    for prot_acc in prot_list:
        bare = prot_acc.split(".")[0]  # strip version for comparison
        found = False

        for record in records:
            for feat in record.features:
                if feat.type != "CDS" or "protein_id" not in feat.qualifiers:
                    continue
                for pid in feat.qualifiers["protein_id"]:
                    if bare == pid.split(".")[0]:
                        hit_start = int(feat.location.start) + 1
                        hit_end   = int(feat.location.end)
                        start = max(1, hit_start - window)
                        end   = min(len(record), hit_end + window)

                        sub = extract_region(record, start, end)
                        rid = safe_record_id(record)
                        outname = os.path.join(
                            output_dir,
                            f"{accession_label}_{rid}_{prot_acc}_{start}_{end}.gbk"
                        )
                        SeqIO.write(sub, outname, "genbank")
                        print(f"  ✔ {prot_acc} → {outname}")
                        found = True
                        break
                if found:
                    break

        if not found:
            print(f"  ✘ Protein accession not found: {prot_acc}")


# ─────────────────────────────────────────────
# MODE C — GENE / PRODUCT NAME SEARCH
# ─────────────────────────────────────────────

SEARCH_QUALIFIERS = ["gene", "locus_tag", "product"]


def match_gene_name(feat, query_lower):
    """
    Return True if any of the searchable qualifiers contain the query string
    (case-insensitive, partial match).
    """
    for qual in SEARCH_QUALIFIERS:
        if qual in feat.qualifiers:
            for val in feat.qualifiers[qual]:
                if query_lower in val.lower():
                    return True
    return False


def find_gbff_files_flat(directory):
    """
    Find all .gbff files directly inside a single flat directory (non-recursive).
    This is the expected layout after downloading with the built-in downloader.
    """
    return sorted(glob.glob(os.path.join(directory, "*.gbff")))


def run_mode_c(output_dir):
    """
    Search a flat directory of .gbff files for a gene/product name,
    extract flanking regions, and write one .gbk per hit.
    Output filenames encode genome + contig + hit index so clinker can tell them apart.

    Expects a flat directory where each file is one genome's .gbff —
    i.e. the directory produced by the built-in NCBI downloader.
    """
    print("\n─── Mode C: Search by Gene / Product Name (multi-genome) ───")
    print("  Expects a flat directory of .gbff files (one per genome).")
    print("  Use the downloader at startup to build this directory if needed.\n")

    gbff_dir = input("  Directory containing .gbff files: ").strip()
    query    = input("  Gene/product name to search (e.g. qseC): ").strip()
    window   = int(input("  Flank size (bp) around each hit: ").strip())
    query_lower = query.lower()

    gbff_files = find_gbff_files_flat(gbff_dir)
    if not gbff_files:
        print(f"\n  [ERROR] No .gbff files found in: {gbff_dir}")
        print("  Make sure the directory contains files ending in .gbff")
        return

    print(f"\n  Found {len(gbff_files)} .gbff file(s). Searching...\n")

    total_hits = 0
    no_hit_genomes = []

    for gbff_path in gbff_files:
        try:
            records = load_all_records(gbff_path)
        except Exception as e:
            print(f"  [WARN] Could not parse {os.path.basename(gbff_path)}: {e}")
            continue

        # Genome label: prefer GCA/GCF accession embedded in filename, else file stem
        stem = os.path.splitext(os.path.basename(gbff_path))[0]
        acc_match = _ACCESSION_RE.search(stem)
        genome_label = acc_match.group(1) if acc_match else stem

        genome_hits = 0

        for record in records:
            contig_label = safe_record_id(record)
            hit_index = 0  # counts hits per contig

            for feat in record.features:
                if feat.type not in ("CDS", "rRNA", "tRNA", "ncRNA"):
                    continue
                if not match_gene_name(feat, query_lower):
                    continue

                hit_index += 1
                hit_start = int(feat.location.start) + 1
                hit_end   = int(feat.location.end)
                start = max(1, hit_start - window)
                end   = min(len(record), hit_end + window)

                sub = extract_region(record, start, end)

                # First hit: no suffix. Additional hits: -2, -3, ...
                suffix = f"-{hit_index}" if hit_index > 1 else ""
                outname = os.path.join(
                    output_dir,
                    f"{genome_label}_{query}{suffix}_{start}_{end}.gbk"
                )
                SeqIO.write(sub, outname, "genbank")
                print(f"  ✔ {genome_label} | contig {contig_label}{suffix} | "
                      f"pos {start}–{end}")
                genome_hits += 1
                total_hits += 1

        if genome_hits == 0:
            no_hit_genomes.append(genome_label)

    # Summary
    print(f"\n  ── Summary ──")
    print(f"  Genomes searched : {len(gbff_files)}")
    print(f"  Genomes with hits: {len(gbff_files) - len(no_hit_genomes)}")
    print(f"  Total hits written: {total_hits}")

    if no_hit_genomes:
        print(f"\n  Genomes with no hits ({len(no_hit_genomes)}):")
        for g in no_hit_genomes:
            print(f"    – {g}")

    if total_hits > 0:
        print(f"\n  Output .gbk files are in: {output_dir}")
        print("  Run clinker with:")
        print(f"    clinker {output_dir}/*.gbk -p {query}_synteny.html\n")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    print("\n╔══════════════════════════════════════════════╗")
    print("║   Genomic Region Extractor + NCBI Downloader  ║")
    print("╚══════════════════════════════════════════════╝\n")

    # ── Optional: download .gbff files first ──
    need_download = input(
        "Do you need to download .gbff files from NCBI first? [y/N]: "
    ).strip().lower()
    if need_download == "y":
        run_downloader()

    # ── Choose extraction mode ──
    print("\nExtraction modes:")
    print("  A — Extract by genomic coordinates (single genome)")
    print("  B — Extract by protein accession(s) (single genome)")
    print("  C — Search by gene/product name across a directory of genomes")
    print("  Q — Quit\n")

    mode = input("Choose mode [A/B/C/Q]: ").strip().upper()

    if mode == "Q":
        print("Bye.")
        return

    if mode in ("A", "B"):
        gbff_path = input("\nPath to .gbff file: ").strip()
        records = load_all_records(gbff_path)
        # Use filename stem as label (e.g. GCF_000123456.1)
        accession_label = os.path.splitext(os.path.basename(gbff_path))[0]
        print(f"Loaded: {accession_label}  ({len(records)} contig(s))")

        output_dir = input("Output directory: ").strip()
        os.makedirs(output_dir, exist_ok=True)

        if mode == "A":
            run_mode_a(records, accession_label, output_dir)
        else:
            run_mode_b(records, accession_label, output_dir)

    elif mode == "C":
        output_dir = input("\nOutput directory for extracted .gbk files: ").strip()
        os.makedirs(output_dir, exist_ok=True)
        run_mode_c(output_dir)

    else:
        print("Invalid option.")
        return

    print("\nAll done.\n")


if __name__ == "__main__":
    main()
