#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import copy
import os

def get_record_by_accession(gbff_path):
    for rec in SeqIO.parse(gbff_path, "genbank"):
        accession = rec.id.split('.')[0]
        return rec, accession
    raise ValueError("No GenBank records found.")

def features_overlap(f_start, f_end, region_start, region_end):
    return (f_end > region_start) and (f_start < region_end)

def shift_location(loc, shift, region_len):
    """Shift and clip FeatureLocation or CompoundLocation."""
    if isinstance(loc, CompoundLocation):
        parts = [shift_location(p, shift, region_len) for p in loc.parts]
        parts = [p for p in parts if p is not None]
        if not parts:
            return None
        return CompoundLocation(parts, loc.operator)

    start = int(loc.start) + shift
    end = int(loc.end) + shift

    if end <= 0 or start >= region_len:
        return None

    if start < 0:
        start = 0
    if end > region_len:
        end = region_len

    return FeatureLocation(start, end, strand=loc.strand)

def extract_region(record, start, end):
    region_start = start - 1
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

        shift = -region_start
        new_loc = shift_location(feat.location, shift, region_len)
        if new_loc is None:
            continue

        new_feat.location = new_loc
        new_features.append(new_feat)

    sub.features = new_features
    return sub

def option_a():
    start = int(input("Enter start coordinate: "))
    end = int(input("Enter end coordinate: "))
    if start > end:
        start, end = end, start
    return start, end

def option_b(record):
    prot_input = input(
        "Enter protein accession(s), comma-separated if multiple: "
    ).strip()
    prot_list = [p.strip() for p in prot_input.split(",")]
    window = int(input("Enter flank size (bp) around each protein: ").strip())

    regions = []
    for prot_acc in prot_list:
        found = False
        for feat in record.features:
            if feat.type == "CDS" and "protein_id" in feat.qualifiers:
                # Check against each item in the list, strip version suffixes
                for pid in feat.qualifiers["protein_id"]:
                    if prot_acc.split('.')[0] == pid.split('.')[0]:
                        hit_start = int(feat.location.start) + 1
                        hit_end = int(feat.location.end)
                        start = max(1, hit_start - window)
                        end = min(len(record), hit_end + window)
                        regions.append((prot_acc, start, end))
                        found = True
                        break
            if found:
                break
        if not found:
            print(f"Warning: Protein accession {prot_acc} not found.")
    return regions

def main():
    print("\n=== Genomic Region Extractor (GeneGraphics-compatible) ===\n")

    gbff_path = input("Path to .gbff file: ").strip()
    record, accession = get_record_by_accession(gbff_path)
    print(f"Loaded genome accession: {accession}")

    output_dir = input("Enter name for output directory: ").strip()
    os.makedirs(output_dir, exist_ok=True)

    mode = input(
        "\nChoose mode:\n"
        "  A = extract by coordinates\n"
        "  B = extract region around protein accession(s)\n"
        "Enter A or B: "
    ).strip().upper()

    if mode == "A":
        start, end = option_a()
        sub = extract_region(record, start, end)
        outname = os.path.join(output_dir, f"{accession}_{start}_{end}.gbk")
        SeqIO.write(sub, outname, "genbank")
        print(f"\n✔ Output written to: {outname}\n")

    elif mode == "B":
        regions = option_b(record)
        for prot_acc, start, end in regions:
            sub = extract_region(record, start, end)
            outname = os.path.join(output_dir, f"{accession}_{prot_acc}_{start}_{end}.gbk")
            SeqIO.write(sub, outname, "genbank")
            print(f"\n✔ Output written to: {outname}")
    else:
        raise ValueError("Invalid option.")

    print("\nAll requested regions processed successfully.\n")

if __name__ == "__main__":
    main()
