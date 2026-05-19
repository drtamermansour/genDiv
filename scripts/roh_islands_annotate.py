"""Identify ROH islands per group and annotate with gene content +
selection-candidate flags.

For each group's per-window ROH frequency landscape produced by
roh_common_landscape.py (landscape.${rg}.tsv):
  1. Compute the top-percentile f_w threshold across all autosomal
     windows (default top 1 %), and additionally flag the absolute
     threshold f_w >= 0.5 ("high-confidence" peak regardless of rank).
  2. Mark "candidate" windows above the percentile threshold.
  3. Merge contiguous candidate windows allowing a gap of up to
     <max-gap-windows> non-candidate windows between them to bridge
     local dips in f_w.
  4. Drop islands narrower than <min-island-kb> kb.
  5. For each surviving island, intersect with the Ensembl EquCab3
     GTF protein-coding gene set and record all overlapping gene
     symbols.
  6. Flag islands that overlap any gene in a curated horse-selection
     candidate list.

A cross-group consolidation pass then bedtools-merges the per-group
island sets into unique regions and records which group(s)
contributed to each region (shared vs gait-specific signals).

When --differential-pair LABEL1:LABEL2 is supplied, the script also
runs a secondary differential scan: per-window delta = f_w(LABEL1)
minus f_w(LABEL2), thresholded at the top-percentile of |delta| (with
positive and negative tails kept separately), then merged + filtered
by minimum island width as above. This identifies regions
specifically enriched in one group relative to the other (e.g., the
DMRT3 gait-keeper locus in Pacer vs Trotter) without requiring those
regions to make the per-group rank-based top-percentile.

Chromosome naming: landscape files use "chrN" / "chrX"; the Ensembl
GTF uses bare "N" / "X". The script strips the "chr" prefix when
matching islands against genes; output retains the landscape naming.

Inputs:
  --group LABEL:LANDSCAPE_TSV     repeat per group
  --gtf                            Ensembl EquCab3 GTF (gzipped OK)
  --candidates                     TSV with at least a 'gene_symbol' column
  --top-percentile                 default 1.0
  --absolute-threshold             default 0.5
  --max-gap-windows                default 2  (=200 kb at 100 kb windows)
  --min-island-kb                  default 500
  --out-dir                        output directory
"""

import argparse
import gzip
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def load_landscape(path):
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str})
    if not {"chrom", "start", "end", "f_w"}.issubset(df.columns):
        raise SystemExit(f"{path} missing required columns (chrom/start/end/f_w)")
    return df.sort_values(["chrom", "start"]).reset_index(drop=True)


def load_candidate_genes(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    if "gene_symbol" not in df.columns:
        raise SystemExit(f"{path} must contain a 'gene_symbol' column")
    return set(df["gene_symbol"].str.strip().str.upper())


def extract_genes_from_gtf(gtf_path):
    """Return DataFrame: chrom, start, end, gene_symbol, gene_biotype."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    name_re = re.compile(r'gene_name "([^"]+)"')
    bio_re = re.compile(r'gene_biotype "([^"]+)"')
    rows = []
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GTF 1-based -> 0-based half-open
            end = int(parts[4])
            attrs = parts[8]
            m_name = name_re.search(attrs)
            m_bio = bio_re.search(attrs)
            rows.append((
                chrom, start, end,
                m_name.group(1) if m_name else "",
                m_bio.group(1) if m_bio else "",
            ))
    return pd.DataFrame(rows, columns=["chrom", "start", "end",
                                       "gene_symbol", "gene_biotype"])


def call_islands(land, top_pct, abs_thresh, max_gap_windows, min_island_kb):
    """Threshold + merge per chromosome. Returns list of island dicts and
    the dynamic top-percentile threshold."""
    thr = float(np.percentile(land["f_w"].to_numpy(), 100.0 - top_pct))
    land = land.copy()
    land["is_candidate"] = land["f_w"] >= thr

    islands = []
    for chrom, grp in land.groupby("chrom", sort=False):
        grp = grp.reset_index(drop=True)
        is_cand = grp["is_candidate"].to_numpy()
        if not is_cand.any():
            continue
        n = len(grp)
        i = 0
        while i < n:
            if not is_cand[i]:
                i += 1
                continue
            start_idx = i
            last_cand_idx = i
            j = i + 1
            while j < n:
                if is_cand[j]:
                    last_cand_idx = j
                    j += 1
                else:
                    # Look ahead up to max_gap_windows non-candidate windows
                    # for the next candidate; if found, absorb the gap.
                    gap_extend = 0
                    k = j
                    while (k < n and not is_cand[k]
                           and gap_extend < max_gap_windows):
                        k += 1
                        gap_extend += 1
                    if k < n and is_cand[k]:
                        j = k
                    else:
                        break

            isl_start = int(grp["start"].iloc[start_idx])
            isl_end = int(grp["end"].iloc[last_cand_idx])
            width_kb = (isl_end - isl_start) / 1000.0
            if width_kb >= min_island_kb:
                island_grp = grp.iloc[start_idx:last_cand_idx + 1]
                peak_row = island_grp.iloc[island_grp["f_w"].argmax()]
                islands.append({
                    "chrom": chrom,
                    "start": isl_start,
                    "end": isl_end,
                    "width_kb": round(width_kb, 1),
                    "n_windows": int(last_cand_idx - start_idx + 1),
                    "peak_f_w": round(float(peak_row["f_w"]), 4),
                    "peak_start": int(peak_row["start"]),
                    "peak_end": int(peak_row["end"]),
                    "high_confidence_abs": float(peak_row["f_w"]) >= abs_thresh,
                })
            i = j if j > start_idx else i + 1

    return islands, thr


def annotate_islands(islands_df, genes_df, candidates):
    """Annotate each island with overlapping protein_coding genes and
    candidate-gene hits. The 'chr' prefix is stripped from island
    chromosomes for the GTF match."""
    if islands_df.empty:
        return islands_df.assign(
            n_genes=pd.Series(dtype=int),
            gene_symbols=pd.Series(dtype=str),
            candidate_genes_hit=pd.Series(dtype=str),
            known_selection_candidate=pd.Series(dtype=str),
        )

    genes_pc = genes_df[genes_df["gene_biotype"] == "protein_coding"].copy()
    genes_pc["chrom_match"] = genes_pc["chrom"].str.replace(r"^chr", "", regex=True)
    isl = islands_df.copy()
    isl["chrom_match"] = isl["chrom"].str.replace(r"^chr", "", regex=True)

    out_rows = []
    for _, row in isl.iterrows():
        hits = genes_pc[
            (genes_pc["chrom_match"] == row["chrom_match"])
            & (genes_pc["end"] > row["start"])
            & (genes_pc["start"] < row["end"])
        ]
        gene_syms = sorted({g for g in hits["gene_symbol"].tolist() if g})
        cand_hits = sorted({g for g in gene_syms if g.upper() in candidates})
        out_rows.append({
            **{k: v for k, v in row.items() if k != "chrom_match"},
            "n_genes": len(gene_syms),
            "gene_symbols": ";".join(gene_syms),
            "candidate_genes_hit": ";".join(cand_hits),
            "known_selection_candidate": "Yes" if cand_hits else "No",
        })
    return pd.DataFrame(out_rows)


def call_differential_islands(land_a, land_b, top_pct, max_gap_windows,
                              min_island_kb):
    """Compute per-window delta = f_w(A) - f_w(B), threshold both tails at
    |delta| >= the top-percentile of |delta|, merge contiguous windows
    allowing gaps, and return (islands_A_enriched, islands_B_enriched,
    threshold).

    Both landscapes must be aligned on (chrom, start, end). Windows
    appearing in only one landscape are dropped from the differential
    set (their delta is undefined for ranking)."""
    merged = land_a[["chrom", "start", "end", "f_w"]].merge(
        land_b[["chrom", "start", "end", "f_w"]],
        on=["chrom", "start", "end"], how="inner",
        suffixes=("_A", "_B"),
    ).sort_values(["chrom", "start"]).reset_index(drop=True)
    merged["delta"] = merged["f_w_A"] - merged["f_w_B"]

    abs_thr = float(np.percentile(merged["delta"].abs().to_numpy(),
                                  100.0 - top_pct))

    def _call_tail(direction):
        """direction = +1 for A-enriched (delta >= +abs_thr),
                       -1 for B-enriched (delta <= -abs_thr)."""
        if direction > 0:
            merged["is_candidate"] = merged["delta"] >= abs_thr
        else:
            merged["is_candidate"] = merged["delta"] <= -abs_thr

        islands = []
        for chrom, grp in merged.groupby("chrom", sort=False):
            grp = grp.reset_index(drop=True)
            is_cand = grp["is_candidate"].to_numpy()
            if not is_cand.any():
                continue
            n = len(grp)
            i = 0
            while i < n:
                if not is_cand[i]:
                    i += 1
                    continue
                start_idx = i
                last_cand_idx = i
                j = i + 1
                while j < n:
                    if is_cand[j]:
                        last_cand_idx = j
                        j += 1
                    else:
                        gap_extend = 0
                        k = j
                        while (k < n and not is_cand[k]
                               and gap_extend < max_gap_windows):
                            k += 1
                            gap_extend += 1
                        if k < n and is_cand[k]:
                            j = k
                        else:
                            break

                isl_start = int(grp["start"].iloc[start_idx])
                isl_end = int(grp["end"].iloc[last_cand_idx])
                width_kb = (isl_end - isl_start) / 1000.0
                if width_kb >= min_island_kb:
                    island_grp = grp.iloc[start_idx:last_cand_idx + 1]
                    if direction > 0:
                        peak_row = island_grp.iloc[
                            island_grp["delta"].argmax()]
                    else:
                        peak_row = island_grp.iloc[
                            island_grp["delta"].argmin()]
                    islands.append({
                        "chrom": chrom,
                        "start": isl_start,
                        "end": isl_end,
                        "width_kb": round(width_kb, 1),
                        "n_windows": int(last_cand_idx - start_idx + 1),
                        "peak_delta": round(float(peak_row["delta"]), 4),
                        "peak_f_w_A": round(float(peak_row["f_w_A"]), 4),
                        "peak_f_w_B": round(float(peak_row["f_w_B"]), 4),
                        "peak_position": int(peak_row["start"]),
                    })
                i = j if j > start_idx else i + 1
        return islands

    return _call_tail(+1), _call_tail(-1), abs_thr


def consolidate_cross_group(per_group_islands, out_path):
    """bedtools-merge all islands across groups, then re-annotate which
    groups contributed to each merged region."""
    all_rows = []
    for label, df in per_group_islands.items():
        if df.empty:
            continue
        for _, r in df.iterrows():
            all_rows.append({"chrom": r["chrom"], "start": int(r["start"]),
                             "end": int(r["end"]), "group": label,
                             "peak_f_w": float(r["peak_f_w"]),
                             "gene_symbols": r["gene_symbols"],
                             "candidate_genes_hit": r["candidate_genes_hit"]})
    if not all_rows:
        print("[islands] no islands found in any group; consolidation skipped")
        return None

    all_df = pd.DataFrame(all_rows)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as tmp:
        all_df[["chrom", "start", "end"]].to_csv(
            tmp.name, sep="\t", index=False, header=False)
        bed_path = tmp.name
    sorted_path = bed_path + ".sorted"
    merged_path = bed_path + ".merged"
    with open(sorted_path, "w") as fh:
        subprocess.run(["sort", "-k1,1", "-k2,2n", bed_path],
                       stdout=fh, check=True)
    with open(merged_path, "w") as fh:
        subprocess.run(["bedtools", "merge", "-i", sorted_path],
                       stdout=fh, check=True)
    merged = pd.read_csv(merged_path, sep="\t", header=None,
                         names=["chrom", "start", "end"], dtype={"chrom": str})
    for p in (bed_path, sorted_path, merged_path):
        Path(p).unlink(missing_ok=True)

    group_labels = sorted(per_group_islands.keys())
    cons_rows = []
    for _, mrow in merged.iterrows():
        region = {
            "chrom": mrow["chrom"],
            "start": int(mrow["start"]),
            "end": int(mrow["end"]),
            "width_kb": round((int(mrow["end"]) - int(mrow["start"])) / 1000.0, 1),
        }
        peaks = {}
        gene_acc, cand_acc = [], []
        for label in group_labels:
            sub = per_group_islands[label]
            if sub.empty:
                continue
            overlap = sub[(sub["chrom"] == region["chrom"])
                          & (sub["end"] > region["start"])
                          & (sub["start"] < region["end"])]
            if not overlap.empty:
                peaks[label] = round(float(overlap["peak_f_w"].max()), 4)
                for gs in overlap["gene_symbols"].tolist():
                    if gs:
                        gene_acc.extend(gs.split(";"))
                for cs in overlap["candidate_genes_hit"].tolist():
                    if cs:
                        cand_acc.extend(cs.split(";"))
        for label in group_labels:
            region[f"present_in_{label}"] = "Yes" if label in peaks else "No"
            region[f"peak_f_w_{label}"] = peaks.get(label, np.nan)
        region["n_groups"] = sum(1 for label in group_labels if label in peaks)
        region["gene_symbols"] = ";".join(sorted(set(g for g in gene_acc if g)))
        region["candidate_genes_hit"] = ";".join(sorted(set(cand_acc)))
        region["known_selection_candidate"] = "Yes" if cand_acc else "No"
        cons_rows.append(region)

    cons_df = pd.DataFrame(cons_rows).sort_values(
        ["chrom", "start"]).reset_index(drop=True)
    cons_df.to_csv(out_path, index=False, float_format="%.4f")
    print(f"[islands] consolidated: {len(cons_df)} unique regions across "
          f"{len(group_labels)} groups; wrote {out_path}")
    return cons_df


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group", action="append", required=True,
                    help="Repeatable: LABEL:LANDSCAPE_TSV.")
    ap.add_argument("--gtf", required=True, help="Ensembl EquCab3 GTF.")
    ap.add_argument("--candidates", required=True,
                    help="TSV with at least a 'gene_symbol' column.")
    ap.add_argument("--top-percentile", type=float, default=1.0)
    ap.add_argument("--absolute-threshold", type=float, default=0.5)
    ap.add_argument("--max-gap-windows", type=int, default=2)
    ap.add_argument("--min-island-kb", type=float, default=500)
    ap.add_argument("--differential-pair", default=None,
                    help="LABEL1:LABEL2 — also run a differential scan on "
                         "f_w(LABEL1) - f_w(LABEL2). Identifies regions "
                         "specifically enriched in one group relative to "
                         "the other.")
    ap.add_argument("--differential-top-percentile", type=float, default=1.0,
                    help="Top-percentile of |delta f_w| to call differential "
                         "candidate windows (default 1.0).")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = load_candidate_genes(args.candidates)
    print(f"[islands] loaded {len(candidates)} candidate gene symbols")
    genes_df = extract_genes_from_gtf(args.gtf)
    pc_count = int((genes_df["gene_biotype"] == "protein_coding").sum())
    print(f"[islands] loaded {len(genes_df)} gene records from GTF "
          f"({pc_count} protein_coding)")

    per_group_islands = {}
    for spec in args.group:
        if ":" not in spec:
            raise SystemExit(f"--group expects LABEL:PATH; got {spec!r}")
        label, path = spec.split(":", 1)
        land = load_landscape(path)
        islands, thr = call_islands(land, args.top_percentile,
                                    args.absolute_threshold,
                                    args.max_gap_windows,
                                    args.min_island_kb)
        islands_df = pd.DataFrame(islands)
        annotated = annotate_islands(islands_df, genes_df, candidates)
        if not annotated.empty:
            annotated = annotated.sort_values(
                ["chrom", "start"]).reset_index(drop=True)
        out_path = out_dir / f"roh_islands.{label}.csv"
        annotated.to_csv(out_path, index=False)
        n_hits = int((annotated["known_selection_candidate"] == "Yes").sum()) \
            if not annotated.empty else 0
        print(f"[islands] {label}: top-{args.top_percentile:g}% f_w threshold "
              f"= {thr:.3f}; {len(annotated)} islands "
              f"(>= {args.min_island_kb} kb); {n_hits} hit known candidates; "
              f"wrote {out_path}")
        per_group_islands[label] = annotated

    consolidate_cross_group(per_group_islands,
                            out_dir / "roh_islands.consolidated.csv")

    # ---------------- differential scan (optional) ----------------
    if args.differential_pair:
        if ":" not in args.differential_pair:
            raise SystemExit("--differential-pair expects LABEL1:LABEL2; "
                             f"got {args.differential_pair!r}")
        a_label, b_label = args.differential_pair.split(":", 1)
        # Reload landscapes — keep raw per-window data, not islands.
        landscapes = {}
        for spec in args.group:
            label, path = spec.split(":", 1)
            if label in (a_label, b_label):
                landscapes[label] = load_landscape(path)
        missing = [l for l in (a_label, b_label) if l not in landscapes]
        if missing:
            raise SystemExit(f"--differential-pair labels {missing} not found "
                             f"in --group inputs")
        a_islands, b_islands, abs_thr = call_differential_islands(
            landscapes[a_label], landscapes[b_label],
            args.differential_top_percentile, args.max_gap_windows,
            args.min_island_kb)
        # Tag the direction column then annotate.
        a_df = pd.DataFrame(a_islands)
        if not a_df.empty:
            a_df.insert(0, "direction", f"{a_label}_enriched")
        b_df = pd.DataFrame(b_islands)
        if not b_df.empty:
            b_df.insert(0, "direction", f"{b_label}_enriched")
        diff_df = pd.concat([a_df, b_df], ignore_index=True)
        if diff_df.empty:
            print(f"[islands] differential {a_label} vs {b_label}: "
                  f"no islands at top-{args.differential_top_percentile:g}% "
                  f"|delta| >= {abs_thr:.3f}")
        else:
            annotated = annotate_islands(diff_df, genes_df, candidates)
            annotated = annotated.sort_values(
                ["direction", "chrom", "start"]).reset_index(drop=True)
            out_path = (out_dir
                        / f"roh_islands.differential.{a_label}_vs_{b_label}.csv")
            annotated.to_csv(out_path, index=False)
            n_a = int((annotated["direction"] == f"{a_label}_enriched").sum())
            n_b = int((annotated["direction"] == f"{b_label}_enriched").sum())
            n_hits = int(
                (annotated["known_selection_candidate"] == "Yes").sum())
            print(f"[islands] differential {a_label} vs {b_label}: "
                  f"top-{args.differential_top_percentile:g}% |delta| threshold "
                  f"= {abs_thr:.3f}; {n_a} {a_label}-enriched + "
                  f"{n_b} {b_label}-enriched islands; "
                  f"{n_hits} hit known candidates; wrote {out_path}")


if __name__ == "__main__":
    main()
