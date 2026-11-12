
#!/usr/bin/env python3

import argparse
import csv
import math
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

    parser.add_argument(
        "--reference-fai",
        help=(
            "FASTA index of the reference assembly; chromosome lengths "
            "for the figure. Required with --plot."
        )
    )

    parser.add_argument(
        "--plot",
        action="store_true",
        help=(
            "Also draw the A3 pairwise ROH comparison figure "
            "(<out>.pairwise_roh_comparison.pdf)."
        )
    )

    parser.add_argument(
        "--plot-png",
        action="store_true",
        help="With --plot, also write a 300 dpi PNG of the figure."
    )

    args = parser.parse_args()

    if args.plot and not args.reference_fai:
        parser.error("--plot requires --reference-fai")

    # Fail before the VCF pass, not after it.
    if args.plot and not os.path.isfile(args.reference_fai):
        parser.error(f"--reference-fai not found: {args.reference_fai}")

    if args.plot_png and not args.plot:
        parser.error("--plot-png requires --plot")

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


def format_segments(segments):
    """Format ROH segments as 'start-end', ';'-joined ('' if none)."""

    return ";".join(f"{start}-{end}" for start, end in segments)


def parse_segments(text):
    """Inverse of format_segments."""

    segments = []

    for item in text.split(";"):
        if item:
            start, end = item.split("-")
            segments.append((int(start), int(end)))

    return segments


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

            # The original ROH segment(s) each animal contributes to this
            # interval, kept so every partition traces back to its parent.
            parents_a = [
                (rstart, rend) for rstart, rend in records_a
                if rstart <= start and rend >= end
            ]

            parents_b = [
                (rstart, rend) for rstart, rend in records_b
                if rstart <= start and rend >= end
            ]

            a_has_roh = bool(parents_a)
            b_has_roh = bool(parents_b)

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
                "roh_a": format_segments(parents_a),
                "roh_b": format_segments(parents_b),
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


# =============================================================
# Optional figure: A3 pairwise ROH comparison
# =============================================================

N_AUTOSOMES = 31

# Original-ROH length classes (Mb), as in summary_nseg_bins.py:
# [1, 3), [3, 5), [5, 10), [10, inf).
ROH_MIN_MB = 1.0
ROH_CLASS_EDGES_MB = (3.0, 5.0, 10.0)
ROH_CLASS_LABELS = ("1 to <3 Mb", "3 to <5 Mb", "5 to <10 Mb", "\u226510 Mb")
ROH_CLASS_COLORS = ("#FECC5C", "#FD8D3C", "#F03B20", "#BD0026")  # YlOrRd

TRACK_BG = "#E5E5E5"
SHARED_COLOR = "#000000"
UNRESOLVED_FACE = "#B3B3B3"
UNRESOLVED_HATCH = "////"
UNRESOLVED_HATCH_COLOR = "#4D4D4D"

BOTH_ROH = "both"
EVALUABLE_STATUSES = {"SHARED_HAPLOTYPE", "DIFFERENT_HAPLOTYPES"}
UNRESOLVED_STATUS = "INSUFFICIENT_JOINT_HOM_SNPs"

# Page geometry in inches (A3 landscape).
PAGE_W, PAGE_H = 16.5, 11.7
MARGIN_L, MARGIN_R = 0.25, 0.35
MARGIN_TOP, MARGIN_BOTTOM = 1.25, 0.35
LABEL_W, COL_GAP = 0.6, 0.35
N_ROWS = 16

# Track y-positions inside a chromosome cell (cell spans 0..1).
TRACK_H = 0.17
TRACK_Y = {"A": 0.635, "B": 0.415, "C": 0.195}


def autosome_number(name):
    """Return 1..31 for an equine autosome name (chr1, ECA1, 1), else None."""

    core = name

    for prefix in ("chr", "eca"):
        if core.lower().startswith(prefix):
            core = core[len(prefix):]
            break

    if core.isdigit() and 1 <= int(core) <= N_AUTOSOMES:
        return int(core)

    return None


def read_autosome_lengths(fai_path):
    """Autosome lengths (bp) keyed 1..31 from a FASTA index."""

    lengths = {}

    with open(fai_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 2:
                raise ValueError(f"Malformed FASTA index line in {fai_path}")

            number = autosome_number(fields[0])

            if number is None:
                continue

            if number in lengths:
                raise ValueError(
                    f"ECA{number} appears twice in {fai_path}"
                )

            lengths[number] = int(fields[1])

    missing = [n for n in range(1, N_AUTOSOMES + 1) if n not in lengths]

    if missing:
        raise ValueError(
            f"{fai_path} lacks autosome(s): "
            + ", ".join(f"ECA{n}" for n in missing)
        )

    return lengths


def autosomal_segments(rohs, label, lengths):
    """
    Group one animal's original ROHs by autosome, checking that each
    lies inside its chromosome and that no two overlap.
    """

    by_chrom = defaultdict(list)
    skipped = 0

    for chrom, start, end in rohs:

        number = autosome_number(chrom)

        if number is None:
            skipped += 1
            continue

        if start < 1 or end > lengths[number]:
            raise ValueError(
                f"{label} ROH {chrom}:{start}-{end} falls outside "
                f"ECA{number} (length {lengths[number]})"
            )

        by_chrom[number].append((start, end))

    for number, segments in by_chrom.items():
        segments.sort()

        for (_, prev_end), (start, _) in zip(segments, segments[1:]):
            if start <= prev_end:
                raise ValueError(
                    f"{label} has overlapping ROHs on ECA{number} "
                    f"near {start}; lengths would be double-counted"
                )

    if skipped:
        print(f"Plot: ignoring {skipped} non-autosomal ROH(s) of {label}")

    return by_chrom


def read_overlap_intervals(tsv_path, match_threshold, lengths,
                           segments_a, segments_b):
    """
    Read the partitioned intervals where both animals are in ROH from
    the pairwise TSV, validating them against the original ROHs.
    Concordance is taken as written, never recomputed.
    """

    metadata = {}
    table_lines = []

    with open(tsv_path) as handle:
        for line in handle:
            if line.startswith("##"):
                key, _, value = line[2:].rstrip("\n").partition("=")
                metadata[key] = value
            else:
                table_lines.append(line)

    tsv_threshold = metadata.get("match_threshold")

    if tsv_threshold is not None and \
            not math.isclose(float(tsv_threshold), match_threshold):
        raise ValueError(
            f"{tsv_path} was classified at match_threshold="
            f"{tsv_threshold}, not {match_threshold}"
        )

    reader = csv.DictReader(table_lines, delimiter="\t")

    required = {
        "chrom", "start", "end", "n_joint_hom", "n_match",
        "concordance", "comparison_status",
    }
    missing = required - set(reader.fieldnames or [])

    if missing:
        raise ValueError(
            f"{tsv_path} lacks column(s): {', '.join(sorted(missing))}"
        )

    has_parents = {"roh_A", "roh_B"} <= set(reader.fieldnames)

    if not has_parents:
        print("Plot: TSV has no roh_A/roh_B columns; checking overlap "
              "intervals against the original ROHs by containment only")

    def check_parent(segments, column, row, number, start, end, label):
        containing = [
            seg for seg in segments.get(number, [])
            if seg[0] <= start and seg[1] >= end
        ]

        if not containing:
            raise ValueError(
                f"{tsv_path} interval {row['chrom']}:{start}-{end} is "
                f"marked as a {label} ROH but lies inside no original "
                f"{label} ROH"
            )

        if has_parents and set(parse_segments(row[column])) != \
                set(containing):
            raise ValueError(
                f"{tsv_path} interval {row['chrom']}:{start}-{end}: "
                f"{column}={row[column]!r} does not match the original "
                f"{label} ROH(s)"
            )

    intervals = []

    for row in reader:

        status = row["comparison_status"]

        if status not in EVALUABLE_STATUSES and status != UNRESOLVED_STATUS:
            continue  # only one animal is in ROH: not drawn

        number = autosome_number(row["chrom"])

        if number is None:
            continue

        start = int(row["start"])
        end = int(row["end"])

        if not 1 <= start <= end <= lengths[number]:
            raise ValueError(
                f"{tsv_path} interval {row['chrom']}:{start}-{end} has "
                f"invalid coordinates for ECA{number}"
            )

        check_parent(segments_a, "roh_A", row, number, start, end, "Horse A")
        check_parent(segments_b, "roh_B", row, number, start, end, "Horse B")

        n_joint_hom = int(row["n_joint_hom"])
        concordance = None
        n_match = None

        if status in EVALUABLE_STATUSES:

            if not row["concordance"] or not row["n_match"]:
                raise ValueError(
                    f"{tsv_path} interval {row['chrom']}:{start}-{end} is "
                    f"{status} but has no concordance"
                )

            concordance = float(row["concordance"])
            n_match = int(row["n_match"])

            if not 0 <= concordance <= 1 or n_joint_hom < 1:
                raise ValueError(
                    f"{tsv_path} interval {row['chrom']}:{start}-{end} has "
                    f"concordance {concordance} from {n_joint_hom} SNPs"
                )

            shared = concordance >= match_threshold

            if shared != (status == "SHARED_HAPLOTYPE"):
                raise ValueError(
                    f"{tsv_path} interval {row['chrom']}:{start}-{end} is "
                    f"{status} at concordance {concordance}, inconsistent "
                    f"with threshold {match_threshold}"
                )

        intervals.append({
            "chrom": number,
            "start": start,
            "end": end,
            "status": status,
            "concordance": concordance,
            "n_joint_hom": n_joint_hom,
            "n_match": n_match,
        })

    return intervals, metadata


def physical_overlap_bp(segments_a, segments_b):
    """
    Base pairs where both animals are in ROH, intersected directly from
    the original segments (sorted and non-overlapping per animal), as an
    independent check on the partitioned intervals.
    """

    total = 0

    for number in set(segments_a) & set(segments_b):

        a = segments_a[number]
        b = segments_b[number]
        i = j = 0

        while i < len(a) and j < len(b):
            start = max(a[i][0], b[j][0])
            end = min(a[i][1], b[j][1])

            if start <= end:
                total += end - start + 1

            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1

    return total


def roh_class_color(start, end):
    """
    Colour of an original ROH by its full length class, or None below
    ROH_MIN_MB (outside every class, so not coloured).
    """

    length_mb = (end - start + 1) / 1e6

    if length_mb < ROH_MIN_MB:
        return None

    for edge, color in zip(ROH_CLASS_EDGES_MB, ROH_CLASS_COLORS):
        if length_mb < edge:
            return color

    return ROH_CLASS_COLORS[-1]


def bar_verts(start, end, y):
    """Rectangle for a 1-based inclusive [start, end] interval, in Mb."""

    x0 = (start - 1) / 1e6
    x1 = end / 1e6

    return [(x0, y), (x1, y), (x1, y + TRACK_H), (x0, y + TRACK_H)]


def summarise_comparison(segments_a, segments_b, overlaps):
    """Summary statistics for the figure title block."""

    def total_mb(segments):
        return sum(
            end - start + 1
            for chrom_segments in segments.values()
            for start, end in chrom_segments
        ) / 1e6

    def interval_mb(statuses):
        return sum(
            iv["end"] - iv["start"] + 1
            for iv in overlaps if iv["status"] in statuses
        ) / 1e6

    evaluable = [iv for iv in overlaps if iv["status"] in EVALUABLE_STATUSES]
    joint = sum(iv["n_joint_hom"] for iv in evaluable)

    return {
        "n_roh_a": sum(len(s) for s in segments_a.values()),
        "n_roh_b": sum(len(s) for s in segments_b.values()),
        "roh_mb_a": total_mb(segments_a),
        "roh_mb_b": total_mb(segments_b),
        "overlap_mb": interval_mb(
            EVALUABLE_STATUSES | {UNRESOLVED_STATUS}
        ),
        "shared_mb": interval_mb({"SHARED_HAPLOTYPE"}),
        "different_mb": interval_mb({"DIFFERENT_HAPLOTYPES"}),
        "n_evaluable": len(evaluable),
        "n_joint_hom": joint,
        # Sum of matches over sum of jointly homozygous SNPs, so each
        # interval counts in proportion to the SNPs behind its estimate.
        "weighted_concordance": (
            sum(iv["n_match"] for iv in evaluable) / joint
            if joint else None
        ),
        "n_unresolved": sum(
            iv["status"] == UNRESOLVED_STATUS for iv in overlaps
        ),
    }


def plot_roh_comparison(tsv_path, rohs_a, rohs_b, animal_a, animal_b,
                        fai_path, match_threshold, out_pdf,
                        out_png=None, png_dpi=300):
    """
    Single-page A3 landscape figure of two animals' ROHs on all 31
    autosomes, drawn to one physical scale, with a third track giving
    the haplotype concordance of their overlapping ROH.

    The horse tracks use the original ROH segments (rohs_a, rohs_b);
    the concordance track uses only the partitioned intervals in
    tsv_path, with concordance exactly as written there.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import ListedColormap, Normalize
    from matplotlib.patches import Rectangle

    lengths = read_autosome_lengths(fai_path)
    segments_a = autosomal_segments(rohs_a, "Horse A", lengths)
    segments_b = autosomal_segments(rohs_b, "Horse B", lengths)
    overlaps, metadata = read_overlap_intervals(
        tsv_path, match_threshold, lengths, segments_a, segments_b
    )
    stats = summarise_comparison(segments_a, segments_b, overlaps)

    partitioned_bp = sum(iv["end"] - iv["start"] + 1 for iv in overlaps)
    direct_bp = physical_overlap_bp(segments_a, segments_b)

    if partitioned_bp != direct_bp:
        raise ValueError(
            f"{tsv_path}: both-in-ROH intervals cover {partitioned_bp} bp "
            f"but the original ROHs overlap by {direct_bp} bp"
        )

    for label, segments in (("Horse A", segments_a), ("Horse B", segments_b)):
        n_short = sum(
            roh_class_color(s, e) is None
            for chrom_segments in segments.values()
            for s, e in chrom_segments
        )
        if n_short:
            print(f"Plot: {n_short} ROH(s) of {label} under "
                  f"{ROH_MIN_MB:g} Mb are outside the length classes "
                  f"and left uncoloured")

    max_mb = max(lengths.values()) / 1e6

    # Blue for different haplotypes, over [0, threshold) only. The pale
    # end is trimmed so low concordance stays visible on the gray track.
    blues = ListedColormap(
        plt.get_cmap("Blues")(np.linspace(0.30, 0.85, 256))
    )
    blue_norm = Normalize(vmin=0.0, vmax=match_threshold)

    col_w = (PAGE_W - MARGIN_L - MARGIN_R - COL_GAP) / 2
    track_w = col_w - LABEL_W
    row_h = (PAGE_H - MARGIN_TOP - MARGIN_BOTTOM) / N_ROWS

    def cell_rect(col, row, x_from=0.0):
        """Figure-fraction rect of a cell's track area from x_from (0..1)."""
        left = MARGIN_L + col * (col_w + COL_GAP) + LABEL_W
        bottom = PAGE_H - MARGIN_TOP - (row + 1) * row_h
        return [
            (left + x_from * track_w) / PAGE_W,
            bottom / PAGE_H,
            (1 - x_from) * track_w / PAGE_W,
            row_h / PAGE_H,
        ]

    def slot(number):
        return (0, number - 1) if number <= 16 else (1, number - 17)

    rc = {"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 10}

    with matplotlib.rc_context(rc):

        fig = plt.figure(figsize=(PAGE_W, PAGE_H))

        for number in range(1, N_AUTOSOMES + 1):

            ax = fig.add_axes(cell_rect(*slot(number)))
            ax.set_xlim(0, max_mb)
            ax.set_ylim(0, 1)
            ax.axis("off")

            ax.text(-0.008, 0.5, f"ECA{number}", transform=ax.transAxes,
                    ha="right", va="center", fontweight="bold")

            ax.add_collection(PolyCollection(
                [bar_verts(1, lengths[number], y) for y in TRACK_Y.values()],
                facecolors=TRACK_BG, edgecolors="none",
            ))

            for track, segments in (("A", segments_a), ("B", segments_b)):
                coloured = [
                    (s, e, roh_class_color(s, e))
                    for s, e in segments.get(number, [])
                    if roh_class_color(s, e) is not None
                ]
                if coloured:
                    ax.add_collection(PolyCollection(
                        [bar_verts(s, e, TRACK_Y[track])
                         for s, e, _ in coloured],
                        facecolors=[color for _, _, color in coloured],
                        edgecolors="none",
                    ))

            chrom_overlaps = [iv for iv in overlaps if iv["chrom"] == number]
            evaluable = [iv for iv in chrom_overlaps
                         if iv["status"] in EVALUABLE_STATUSES]
            unresolved = [iv for iv in chrom_overlaps
                          if iv["status"] == UNRESOLVED_STATUS]

            if evaluable:
                ax.add_collection(PolyCollection(
                    [bar_verts(iv["start"], iv["end"], TRACK_Y["C"])
                     for iv in evaluable],
                    facecolors=[
                        SHARED_COLOR if iv["status"] == "SHARED_HAPLOTYPE"
                        else blues(blue_norm(iv["concordance"]))
                        for iv in evaluable
                    ],
                    edgecolors="none",
                ))

            if unresolved:
                # Hatch lines take the edge colour; linewidth 0 keeps the
                # outline from widening the segment.
                ax.add_collection(PolyCollection(
                    [bar_verts(iv["start"], iv["end"], TRACK_Y["C"])
                     for iv in unresolved],
                    facecolors=UNRESOLVED_FACE,
                    edgecolors=UNRESOLVED_HATCH_COLOR,
                    linewidths=0, hatch=UNRESOLVED_HATCH,
                ))

        # ---- Coordinate reference: the empty bottom-right slot, drawn
        # on the same axis range as every chromosome cell.
        ruler = fig.add_axes(cell_rect(1, N_ROWS - 1))
        ruler.set_xlim(0, max_mb)
        ruler.set_ylim(0, 1)
        ruler.axis("off")
        ruler.plot([0, max_mb], [0.62, 0.62], color="black", lw=0.8)

        for tick in np.arange(0, max_mb + 1e-9, 50):
            ruler.plot([tick, tick], [0.62, 0.78], color="black", lw=0.8)
            ruler.text(tick, 0.36, f"{tick:.0f}", ha="center", va="center",
                       fontsize=9)

        ruler.text(max_mb, 0.36, "Mb", ha="right", va="center", fontsize=9)
        ruler.text(0, 0.05, "Physical position, common scale for all "
                   "chromosomes", fontsize=9, va="center", color="#404040")

        # ---- Legend: right of the short chromosomes at the foot of
        # column 2, which never reach that far along the shared axis.
        legend_rows = 0
        reach = 0.0

        for row in range(N_ROWS - 2, -1, -1):
            fraction = lengths[row + 17] / 1e6 / max_mb
            if fraction > 0.32:
                break
            legend_rows += 1
            reach = max(reach, fraction)

        if legend_rows < 6:
            raise ValueError("Too few short chromosomes to host the legend")

        x_from = reach + 0.06
        lg_left_in = MARGIN_L + col_w + COL_GAP + LABEL_W + x_from * track_w
        lg_w_in = (1 - x_from) * track_w
        lg_top_in = PAGE_H - MARGIN_TOP - (N_ROWS - 1 - legend_rows) * row_h
        lg_h_in = legend_rows * row_h

        legend = fig.add_axes([
            lg_left_in / PAGE_W, (lg_top_in - lg_h_in) / PAGE_H,
            lg_w_in / PAGE_W, lg_h_in / PAGE_H,
        ])
        legend.set_xlim(0, lg_w_in)
        legend.set_ylim(lg_h_in, 0)  # inches from the top
        legend.axis("off")

        sw_w, sw_h = 0.42, 0.13
        y = 0.15

        def heading(text):
            nonlocal y
            legend.text(0, y, text, fontweight="bold", va="center")
            y += 0.28

        def swatch(x, label, **patch):
            legend.add_patch(Rectangle(
                (x, y - sw_h / 2), sw_w, sw_h, **patch
            ))
            legend.text(x + sw_w + 0.1, y, label, va="center", fontsize=9.5)

        heading("Tracks per chromosome")
        for track, label in (
            ("A", f"Top: Horse A ({animal_a})"),
            ("B", f"Middle: Horse B ({animal_b})"),
            ("C", "Bottom: haplotype concordance"),
        ):
            swatch(0.1, label, facecolor=TRACK_BG, edgecolor="none")
            y += 0.22
        y += 0.2

        heading("ROH length class (full original segment)")
        for i, (label, color) in enumerate(
                zip(ROH_CLASS_LABELS, ROH_CLASS_COLORS)):
            swatch(0.1 + (i % 2) * 2.0, label,
                   facecolor=color, edgecolor="none")
            if i % 2:
                y += 0.24
        y += 0.2

        heading("Concordance where both horses are in ROH")
        legend.text(0.1, y, "Jointly homozygous, non-missing SNPs only",
                    fontsize=9, color="#404040", va="center")
        y += 0.3
        swatch(0.1, f"Shared haplotype (\u2265{match_threshold:.0%})",
               facecolor=SHARED_COLOR, edgecolor="none")
        y += 0.3

        legend.text(0.1, y, f"Different haplotypes (<{match_threshold:.0%})",
                    fontsize=9.5, va="center")
        y += 0.2
        cbar_w_in = min(3.2, lg_w_in - 0.4)
        cax = fig.add_axes([
            (lg_left_in + 0.1) / PAGE_W,
            (lg_top_in - y - sw_h) / PAGE_H,
            cbar_w_in / PAGE_W, sw_h / PAGE_H,
        ])
        cbar = fig.colorbar(
            ScalarMappable(norm=blue_norm, cmap=blues),
            cax=cax, orientation="horizontal",
        )
        ticks = [t for t in (0, 0.25, 0.5, 0.75) if t < match_threshold]
        cbar.set_ticks(ticks + [match_threshold])
        cbar.set_ticklabels([f"{t:.0%}" for t in ticks + [match_threshold]])
        cbar.ax.tick_params(labelsize=8.5, length=2)
        cbar.outline.set_linewidth(0.5)
        y += sw_h + 0.5

        min_snps = metadata.get("min_joint_hom_snps")
        swatch(0.1, (f"Unresolved (<{min_snps} jointly homozygous SNPs)"
                     if min_snps else
                     "Unresolved (too few jointly homozygous SNPs)"),
               facecolor=UNRESOLVED_FACE, edgecolor=UNRESOLVED_HATCH_COLOR,
               linewidth=0, hatch=UNRESOLVED_HATCH)
        y += 0.3
        legend.text(0.1, y, "Only one horse in ROH: no concordance segment",
                    fontsize=9, color="#404040", va="center")

        # ---- Title and summary statistics
        fig.text(0.5, 1 - 0.32 / PAGE_H,
                 f"Pairwise ROH comparison: {animal_a} vs {animal_b}",
                 ha="center", va="center", fontsize=16, fontweight="bold")

        wc = stats["weighted_concordance"]
        line1 = (
            f"Horse A {animal_a}: {stats['n_roh_a']} ROHs, "
            f"{stats['roh_mb_a']:.1f} Mb      "
            f"Horse B {animal_b}: {stats['n_roh_b']} ROHs, "
            f"{stats['roh_mb_b']:.1f} Mb      "
            f"Physical ROH overlap: {stats['overlap_mb']:.1f} Mb"
        )
        line2 = (
            "Weighted concordance: "
            + (f"{wc:.1%} ({stats['n_evaluable']} evaluable intervals, "
               f"{stats['n_joint_hom']} jointly homozygous SNPs)"
               if wc is not None else "n/a (no evaluable intervals)")
            + f"      Shared {stats['shared_mb']:.1f} Mb, "
            f"different {stats['different_mb']:.1f} Mb      "
            f"Unresolved intervals: {stats['n_unresolved']}"
        )
        fig.text(0.5, 1 - 0.68 / PAGE_H, line1, ha="center", va="center",
                 fontsize=11)
        fig.text(0.5, 1 - 0.95 / PAGE_H, line2, ha="center", va="center",
                 fontsize=11)

        fig.savefig(out_pdf)

        if out_png:
            fig.savefig(out_png, dpi=png_dpi)

        plt.close(fig)

    return stats


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
        "roh_A",
        "roh_B",
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
                "roh_A": interval["roh_a"],
                "roh_B": interval["roh_b"],
            }

            writer.writerow(row)

    print("\nComparison complete.")
    print(f"Matching threshold: {args.match_threshold:.2%}")
    print(f"Minimum jointly homozygous SNPs: "
          f"{args.min_joint_hom_snps}")
    print(f"Intervals: {tsv_path}")
    print(f"FASTA: {fasta_path}")

    # ---------------------------------------------------------
    # 6. Optional figure
    # ---------------------------------------------------------

    if args.plot:

        pdf_path = args.out + ".pairwise_roh_comparison.pdf"
        png_path = (
            args.out + ".pairwise_roh_comparison.png"
            if args.plot_png else None
        )

        plot_roh_comparison(
            tsv_path=tsv_path,
            rohs_a=rohs_a,
            rohs_b=rohs_b,
            animal_a=animal_a,
            animal_b=animal_b,
            fai_path=args.reference_fai,
            match_threshold=args.match_threshold,
            out_pdf=pdf_path,
            out_png=png_path,
        )

        print(f"Figure: {pdf_path}")

        if png_path:
            print(f"Figure: {png_path}")


if __name__ == "__main__":
    main()

