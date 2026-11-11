
#!/usr/bin/env python3

import argparse
import csv
import os
from collections import defaultdict
from bisect import bisect_right

from cyvcf2 import VCF


BASES = {"A", "C", "G", "T"}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare ROH haplotypes between two animals using "
            "jointly homozygous SNPs."
        )
    )

    parser.add_argument("--vcf", required=True,
                        help="Common phased VCF (bgzipped and indexed).")

    parser.add_argument("--roh-trotter", required=True,
                        help="BCFtools ROH Trotter output.")

    parser.add_argument("--roh-pacer", required=True,
                        help="BCFtools ROH Pacer output.")

    parser.add_argument("--animal-a", required=True,
                        help="Sample ID of animal A.")

    parser.add_argument("--animal-b", required=True,
                        help="Sample ID of animal B.")

    parser.add_argument("--out", required=True,
                        help="Output prefix.")

    parser.add_argument(
        "--match-threshold",
        type=float,
        default=0.90,
        help=(
            "Minimum genotype concordance for shared haplotype "
            "(default: 0.90)."
        )
    )

    parser.add_argument(
        "--min-joint-hom-snps",
        type=int,
        default=10,
        help=(
            "Minimum number of jointly homozygous SNPs required "
            "for classification (default: 10)."
        )
    )

    args = parser.parse_args()

    if not 0 < args.match_threshold <= 1:
        parser.error("--match-threshold must be greater than 0 and <= 1")

    if args.min_joint_hom_snps < 1:
        parser.error("--min-joint-hom-snps must be at least 1")

    if args.animal_a == args.animal_b:
        parser.error("Animal A and Animal B must be different samples.")

    return args


def read_rohs(roh_files, sample):
    """
    Read RG records for one sample from the supplied BCFtools ROH files.

    Expected RG fields:
    RG Sample Chromosome Start End Length ...
    """

    records = []

    for path in roh_files:
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue

                fields = line.strip().split()

                if len(fields) < 5 or fields[0] != "RG":
                    continue

                if fields[1] != sample:
                    continue

                chrom = fields[2]
                start = int(fields[3])
                end = int(fields[4])

                if start > end:
                    raise ValueError(
                        f"Invalid ROH coordinates in {path}: "
                        f"{chrom}:{start}-{end}"
                    )

                records.append((chrom, start, end))

    return records


def partition_rohs(rohs_a, rohs_b):
    """
    Partition the union of both animals' ROHs into non-overlapping
    intervals defined by all ROH start and end+1 coordinates.

    Each resulting interval has a constant ROH-carrier status.
    """

    by_chrom = defaultdict(lambda: {"A": [], "B": []})

    for chrom, start, end in rohs_a:
        by_chrom[chrom]["A"].append((start, end))

    for chrom, start, end in rohs_b:
        by_chrom[chrom]["B"].append((start, end))

    intervals = []
    interval_number = 0

    for chrom in sorted(by_chrom):

        records_a = by_chrom[chrom]["A"]
        records_b = by_chrom[chrom]["B"]

        boundaries = set()

        for start, end in records_a + records_b:
            boundaries.add(start)
            boundaries.add(end + 1)

        boundaries = sorted(boundaries)

        for i in range(len(boundaries) - 1):

            start = boundaries[i]
            end = boundaries[i + 1] - 1

            if start > end:
                continue

            a_has_roh = any(
                rstart <= start and rend >= end
                for rstart, rend in records_a
            )

            b_has_roh = any(
                rstart <= start and rend >= end
                for rstart, rend in records_b
            )

            # Retain only intervals covered by at least one animal.
            if not a_has_roh and not b_has_roh:
                continue

            interval_number += 1

            intervals.append({
                "interval_id": f"interval_{interval_number:06d}",
                "chrom": chrom,
                "start": start,
                "end": end,
                "a_has_roh": a_has_roh,
                "b_has_roh": b_has_roh,
            })

    return intervals


def classify_genotype(genotype):
    """
    Return genotype state and allele index.

    States:
      HOM      non-missing homozygous genotype
      HET      heterozygous genotype
      MISSING  missing genotype
    """

    allele1 = int(genotype[0])
    allele2 = int(genotype[1])

    if allele1 < 0 or allele2 < 0:
        return "MISSING", None

    if allele1 != allele2:
        return "HET", None

    return "HOM", allele1


def write_fasta_record(handle, interval_id, sample, label, sequence):
    if not sequence:
        return

    header = (
        f">{interval_id}|{sample}|{label}|"
        f"n_homozygous_sites={len(sequence)}"
    )

    handle.write(header + "\n")

    sequence = "".join(sequence)

    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


def main():

    args = parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)),
                exist_ok=True)

    animal_a = args.animal_a
    animal_b = args.animal_b

    # ---------------------------------------------------------
    # 1. Read ROHs for the two animals
    # ---------------------------------------------------------

    roh_files = [args.roh_trotter, args.roh_pacer]

    rohs_a = read_rohs(roh_files, animal_a)
    rohs_b = read_rohs(roh_files, animal_b)

    if not rohs_a:
        raise ValueError(f"No ROHs found for animal A: {animal_a}")

    if not rohs_b:
        raise ValueError(f"No ROHs found for animal B: {animal_b}")

    print(f"{animal_a}: {len(rohs_a)} ROHs")
    print(f"{animal_b}: {len(rohs_b)} ROHs")

    # ---------------------------------------------------------
    # 2. Partition the union of ROH intervals
    # ---------------------------------------------------------

    intervals = partition_rohs(rohs_a, rohs_b)

    print(f"Partitioned intervals: {len(intervals)}")

    # Prepare run metadata for output TSV
    metadata = [
        f"##animal_A={animal_a}",
        f"##animal_B={animal_b}",
        f"##n_ROH_animal_A={len(rohs_a)}",
        f"##n_ROH_animal_B={len(rohs_b)}",
        f"##match_threshold={args.match_threshold}",
        f"##min_joint_hom_snps={args.min_joint_hom_snps}",
    ]

    # Prepare chromosome-specific interval lookup.
    intervals_by_chrom = defaultdict(list)

    for interval in intervals:
        intervals_by_chrom[interval["chrom"]].append(interval)

    interval_starts = {}

    for chrom, chrom_intervals in intervals_by_chrom.items():
        chrom_intervals.sort(key=lambda x: x["start"])
        interval_starts[chrom] = [
            x["start"] for x in chrom_intervals
        ]

    # Initialize statistics and sequence storage.
    for interval in intervals:

        interval.update({
            "n_variant_sites": 0,
            "n_hom_a": 0,
            "n_hom_b": 0,
            "n_het_a": 0,
            "n_het_b": 0,
            "n_missing_a": 0,
            "n_missing_b": 0,
            "n_joint_hom": 0,
            "n_match": 0,
            "n_discordant": 0,
            "seq_a": [],
            "seq_b": [],
        })

    # ---------------------------------------------------------
    # 3. Open VCF and identify samples
    # ---------------------------------------------------------

    vcf = VCF(args.vcf)

    samples = vcf.samples

    if animal_a not in samples:
        raise ValueError(f"Animal A not found in VCF: {animal_a}")

    if animal_b not in samples:
        raise ValueError(f"Animal B not found in VCF: {animal_b}")

    idx_a = samples.index(animal_a)
    idx_b = samples.index(animal_b)

    print("Extracting genotypes from VCF...")

    # ---------------------------------------------------------
    # 4. Extract genotypes for each partitioned interval
    # ---------------------------------------------------------

    for variant in vcf:

        chrom = variant.CHROM
        pos = variant.POS

        if chrom not in intervals_by_chrom:
            continue

        # Retain biallelic SNPs only.
        if variant.REF not in BASES:
            continue

        if len(variant.ALT) != 1:
            continue

        alt = variant.ALT[0]

        if alt not in BASES:
            continue

        chrom_intervals = intervals_by_chrom[chrom]
        starts = interval_starts[chrom]

        interval_index = bisect_right(starts, pos) - 1

        if interval_index < 0:
            continue

        interval = chrom_intervals[interval_index]

        if pos > interval["end"]:
            continue

        a_has_roh = interval["a_has_roh"]
        b_has_roh = interval["b_has_roh"]

        interval["n_variant_sites"] += 1

        gt_a = variant.genotypes[idx_a]
        gt_b = variant.genotypes[idx_b]

        state_a, allele_a = classify_genotype(gt_a)
        state_b, allele_b = classify_genotype(gt_b)

        # Count genotype states only for animals carrying an ROH
        # in this interval.
        if a_has_roh:
            if state_a == "HOM":
                interval["n_hom_a"] += 1
            elif state_a == "HET":
                interval["n_het_a"] += 1
            else:
                interval["n_missing_a"] += 1

        if b_has_roh:
            if state_b == "HOM":
                interval["n_hom_b"] += 1
            elif state_b == "HET":
                interval["n_het_b"] += 1
            else:
                interval["n_missing_b"] += 1

        # Both animals have an ROH:
        # retain only jointly homozygous calls.
        if a_has_roh and b_has_roh:

            if state_a != "HOM" or state_b != "HOM":
                continue

            interval["n_joint_hom"] += 1

            base_a = (
                variant.REF if allele_a == 0 else alt
            )

            base_b = (
                variant.REF if allele_b == 0 else alt
            )

            interval["seq_a"].append(base_a)
            interval["seq_b"].append(base_b)

            if allele_a == allele_b:
                interval["n_match"] += 1
            else:
                interval["n_discordant"] += 1

        # Only animal A has an ROH:
        # extract its non-missing homozygous calls.
        elif a_has_roh:

            if state_a != "HOM":
                continue

            base_a = (
                variant.REF if allele_a == 0 else alt
            )

            interval["seq_a"].append(base_a)

        # Only animal B has an ROH.
        elif b_has_roh:

            if state_b != "HOM":
                continue

            base_b = (
                variant.REF if allele_b == 0 else alt
            )

            interval["seq_b"].append(base_b)

    vcf.close()

    # ---------------------------------------------------------
    # 5. Classify haplotypes and calculate concordance
    # ---------------------------------------------------------

    tsv_path = args.out + ".pairwise_roh_intervals.tsv"
    fasta_path = args.out + ".pairwise_roh_haplotypes.fasta"

    output_columns = [
        "interval_id",
        "chrom",
        "start",
        "end",
        "n_snps",
        "n_het_A",
        "n_het_B",
        "n_missing_A",
        "n_missing_B",
        "n_joint_hom",
        "n_match",
        "concordance",
        "animal_A_hap",
        "animal_B_hap",
        "comparison_status",
    ]


    with open(tsv_path, "w", newline="") as tsv_handle, \
            open(fasta_path, "w") as fasta_handle:

        # Metadata header lines
        for line in metadata:
            tsv_handle.write(line + "\n")

        writer = csv.DictWriter(
            tsv_handle,
            fieldnames=output_columns,
            delimiter="\t"
        )

        writer.writeheader()

        for interval_number, interval in enumerate(intervals, start=1):

            a_has_roh = interval["a_has_roh"]
            b_has_roh = interval["b_has_roh"]

            n_compared = 0
            concordance = None

            hap_a = "No_ROH"
            hap_b = "No_ROH"

            # Both animals carry an ROH.
            if a_has_roh and b_has_roh:

                n_compared = interval["n_joint_hom"]

                if n_compared < args.min_joint_hom_snps:

                    hap_a = "UNRESOLVED"
                    hap_b = "UNRESOLVED"

                    status = "INSUFFICIENT_JOINT_HOM_SNPs"

                else:

                    n_match = interval["n_match"]

                    concordance = n_match / n_compared

                    if concordance >= args.match_threshold:

                        hap_a = "HAP1"
                        hap_b = "HAP1"

                        status = "SHARED_HAPLOTYPE"

                    else:

                        hap_a = "HAP1"
                        hap_b = "HAP2"

                        status = "DIFFERENT_HAPLOTYPES"

            # Only animal A carries an ROH.
            elif a_has_roh:

                n_compared = interval["n_hom_a"]

                if n_compared >= args.min_joint_hom_snps:
                    hap_a = "HAP1"
                else:
                    hap_a = "UNRESOLVED"

                status = "ONLY_A_ROH"

            # Only animal B carries an ROH.
            else:

                n_compared = interval["n_hom_b"]

                if n_compared >= args.min_joint_hom_snps:
                    hap_b = "HAP1"
                else:
                    hap_b = "UNRESOLVED"

                status = "ONLY_B_ROH"

            # FASTA output
            write_fasta_record(
                fasta_handle,
                str(interval_number),
                animal_a,
                hap_a,
                interval["seq_a"]
            )

            write_fasta_record(
                fasta_handle,
                str(interval_number),
                animal_b,
                hap_b,
                interval["seq_b"]
            )

            row = {
                "interval_id": interval_number,
                "chrom": interval["chrom"],
                "start": interval["start"],
                "end": interval["end"],
                "n_snps": interval["n_variant_sites"],
                "n_het_A": interval["n_het_a"],
                "n_het_B": interval["n_het_b"],
                "n_missing_A": interval["n_missing_a"],
                "n_missing_B": interval["n_missing_b"],
                "n_joint_hom": interval["n_joint_hom"],
                "n_match": (
                    interval["n_match"]
                    if a_has_roh and b_has_roh else ""
                ),
                "concordance": (
                    f"{concordance:.6f}"
                    if concordance is not None else ""
                ),
                "animal_A_hap": hap_a,
                "animal_B_hap": hap_b,
                "comparison_status": status,
            }

            writer.writerow(row)

    print("\nComparison complete.")
    print(f"Matching threshold: {args.match_threshold:.2%}")
    print(f"Minimum jointly homozygous SNPs: "
          f"{args.min_joint_hom_snps}")
    print(f"Intervals: {tsv_path}")
    print(f"FASTA: {fasta_path}")


if __name__ == "__main__":
    main()

