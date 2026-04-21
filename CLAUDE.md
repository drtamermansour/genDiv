# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a genomic diversity assessment pipeline for Standardbred horses (~500+ individuals). It analyzes SNP genotyping data to compute genetic diversity metrics (ROH, heterozygosity, relatedness, PCA) and generate population-level reference ranges for breeding programs.

See also:
- `README.md` — human-oriented setup recipe (conda/mamba environment).
- `MIGRATION.md` — the filename / schema contract between this pipeline and the downstream GPA report pipeline (`../GPA/`). Update it whenever a producer here or a consumer there changes.
- `scripts/validate_popRefs.sh` — dual-mode validator for the per-group reference files; runs at end of the pipeline (and at end of GPA's `create_popFiles.sh` once the GPA-side PR lands).
- `scripts/benchmark/run_benchmark.sh` — builds a 25-Trotter + 25-Pacer subset from an existing full-run `OUTPUT_DIR` and runs `per_group.sh × 3 + aggregate.sh` end-to-end for fast regression testing.

## Running the Pipeline

The pipeline is a thin wrapper that composes three subscripts:

```bash
conda activate genDiv
bash ./genDiversity.sh
```

The wrapper sources shared CONFIG/helpers, invokes shared preprocessing once, loops the per-group stage three times (wholePop/Trotter/Pacer), then runs cross-group aggregation. Each subscript is independently runnable, which is handy for debugging or refreshing one stage:

```bash
bash ./genDiversity_shared.sh                 # preprocessing + whole-pop KING / FST / A_e
bash ./genDiversity_per_group.sh wholePop     # per-group metrics (repeat for Trotter, Pacer)
bash ./genDiversity_per_group.sh Trotter
bash ./genDiversity_per_group.sh Pacer
bash ./genDiversity_aggregate.sh              # cross-group summaries + plots
```

## Environment Setup

Uses conda/mamba with a named environment `genDiv`. Full setup command lives in `README.md`; below is a summary focused on what the pipeline actually invokes at runtime.

**Runtime dependencies called by `genDiversity.sh`:**
- CLI tools: `plink`, `plink2`, `bcftools`, `beagle`, `rclone`
- Interpreters: `python` (with `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`), `Rscript` (with `ggplot2`, `gridExtra`, `viridis`, `reshape2`, `GGally`, `effsize`)

**In the env recipe but not currently invoked by the pipeline:** `gcta`, `bedtools`, `snakemake`.

## Architecture

### Entry Point

The pipeline is split across five bash files at the repo root:

| File | Role | Invocation |
|---|---|---|
| `genDiversity.sh` | Thin wrapper: sources common, runs shared → loops per_group × 3 → runs aggregate | `bash genDiversity.sh` |
| `genDiversity_common.sh` | CONFIG, helpers (`log`, `run_python`, `run_r`, `upload`), `ERR` trap, log redirection (guarded by `GENDIV_LOG_SETUP` so subscripts don't double-log) | sourced |
| `genDiversity_shared.sh` | Whole-pop preprocessing (Sections 1–3) + A_e + FST + effective-genome-length + autosomes.genome + whole-pop KING + IBS + `samples.{wholePop,Trotter,Pacer}.txt` + `sample_groups.tsv` | `bash genDiversity_shared.sh` |
| `genDiversity_per_group.sh` | Per-group metric stage; takes `$rg ∈ {wholePop, Trotter, Pacer}`, runs 3× | `bash genDiversity_per_group.sh <rg>` |
| `genDiversity_aggregate.sh` | Cross-group summaries (twoGait / threeBooksize), F_ROH histograms + gait/book-size stratified summaries, Froh-vs-ROHsh plots, merged-kin top-pair cross-reference | `bash genDiversity_aggregate.sh` |

All scripts respect `OUTPUT_DIR` as an env override, so re-running into an existing folder (e.g., to refresh just one group) works the same way as a fresh timestamped run.

`set -eo pipefail` applies everywhere; the `ERR` trap reports the failing line number from whichever subscript errored.

Conceptual workflow phases:

1. **Data Download & Preprocessing** (`shared.sh`) — rclone from Google Drive, PLINK ID updates, EquCab3 remapping, SNP deduplication, PLINK1→PLINK2 conversion.
2. **Data Exploration** (`shared.sh`) — sex validation (X chr F-stats), PAR removal, HWE analysis.
3. **Final Filtering** (`shared.sh`) — apply missingness/MAF/HWE thresholds to produce the clean dataset; LD pruning.
4. **Whole-pop diversity metrics** (`shared.sh`) — A_e (effective-allele number, per SNP and stratified by gait and gait×book-size via `--loop-cats`), FST between subpopulations (sex / gait / book-size), effective autosomal genome length, autosomes.genome.
5. **Whole-pop KING-robust kinship + IBS** (`shared.sh`) — plink2 `--make-king-table`, IBS augmentation of `.kin0`, first-degree-pair filter (`related`), KING-vs-IBS correlation plot.
6. **Per-group reference files** (`per_group.sh`, 3×) — 15 sections mirroring the original pipeline order. Highlights: per-group afreq, PCA + overlays (wSex / wGait wholePop-only, wBook_Size / wCOI all groups), FST book-size-within-gait, per-group F_SNP `.het`, per-group bcftools roh + L1/L2/L3 + consensus ROH (nested over book-size for Trotter/Pacer), F_ROH summary, per-group GRM + ROHRM + analysis_comparison across all `$ROH_CUTOFFS`, `related_${rg}` filter, PCA pairwise Euclidean + KING merge + correlation plots, cross-method correlation plot, and COI-vs-F_ROH / F_ROH-vs-D_ROH / F_SNP-vs-D_ROH / F_SNP-vs-F_ROH doubleAnn plots.
7. **Cross-group aggregation** (`aggregate.sh`) — twoGait and threeBooksize per-sample ROH_sh concatenations, F_ROH histograms + `roh_high.csv` + gait / book-size stratified F_ROH summaries, Froh-vs-ROHsh plots iterating wholePop / twoGait / threeBooksize, merged-kin-sorted-top cross-reference.
8. **Upload** — `rclone` is invoked throughout each subscript; there's no single upload phase.

### Python Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `ROHRM_Creator.py` | Builds the ROH-based Relationship Matrix (ROHRM) — implements Howard et al. C++ logic in Python using window-based haplotype matching on phased VCF data |
| `analysis_comparison.py` | Compares ROHRM vs standard GRM; `RobustMatrixComparator` class |
| `summary_roh.py` | ROH summary statistics grouped by subpopulation. `-n/--n-numeric-cols N` (default 3) controls how many numeric columns are summarised; the column right after the numeric block is the grouping factor |
| `summary_het.py` | Heterozygosity summary stats (observed/expected, F-coefficients) |
| `roh_plot.py` | Scatter plots correlating F_ROH vs consensus ROH sharing |
| `roh_histograms.py` | Unified histogram script: use `--metric ratio` (ROH_shared/F_ROH) or `--metric shared` (Percent_of_Consensus_ROH) |
| `utils.py` | Shared constants (`BOOK_SIZE_ORDER`, `BOOK_SIZE_COLORS`) and `format_stats()` used by summary scripts |

All Python scripts use `argparse`; run with `--help` to see usage.

### R Scripts (`scripts/`)

Grouped by role. Each script's `genDiversity.sh` call site in parentheses.

**PCA:**

| Script | Purpose |
|--------|---------|
| `pca_plots.R` | Unified PCA plots: positional args `<prefix> <eigenvec_suffix> <color_col> <output.png> [n_pcs=3] [color_type=factor\|numeric]` |

**Correlation heatmaps & scatterplot matrices** — one script, four modes (formerly `correlation_plot_multiway{,_v2,_v2e,_v3}.R`):

| Script + mode | Inputs | Purpose |
|---------------|--------|---------|
| `correlation_plot.R --mode basic` | `<homo_file> <het_file>` | KB, KBAVG, HET_diff, F (individual-level) |
| `correlation_plot.R --mode froh` | `<diag_file> <het_file> <froh_file>` | F_SNP, F_ROH, D_ROH, D_STD |
| `correlation_plot.R --mode froh-cons` | `<diag_file> <het_file> <froh_file> <cons_file>` | Adds ROH_sh (consensus share) |
| `correlation_plot.R --mode pairwise` | `<diag_file> <kingkin_file> <eucl_file>` | Pair-level: G_STD, G_ROH, King_kin, IBS |

**Scatter plots:**

| Script | Purpose |
|--------|---------|
| `plot_correlation.R` | Scatter of two columns with single-column color |
| `plot_correlation_withColorsAndShapes.R` | Scatter with both color and shape mappings |

**Statistics / plotting scripts** (hard-coded input paths rooted at `results_<timestamp>/`, no CLI args):

| Script | Purpose |
|--------|---------|
| `fst_stats.R` | FST statistics across subpopulations; stdout captured to `results_<timestamp>/divStats/fst_stats.txt`; writes `results_<timestamp>/divStats/Fst_Analysis_Results_Adjusted.csv` |
| `plot_Ae.R` | Violin plots of effective allele number (Ae) by gait and book size |
| `effAllele_stats.R` | Effect size analysis for Ae; stdout captured to `results_<timestamp>/divStats/effAllele_stats.txt` |

### Data Flow

```
Google Drive (rclone download)
  ↓
genDiversity_shared.sh  (runs once)
  PLINK preprocessing → QC & filtering → LD pruning
  → A_e (whole-pop + gait + gait×book-size) + FST (sex / gait / book-size)
  → effective_autosomal_genome_length.txt + autosomes.genome
  → whole-pop KING + IBS + related + KING-vs-IBS correlation plot
  → preprocess/samples.{wholePop,Trotter,Pacer}.txt + sample_groups.tsv
  ↓
genDiversity_per_group.sh  (runs 3× — wholePop, Trotter, Pacer)
  §1  pruned.${rg}.afreq
  §2  plink2 --pca (+ wSex/wGait/wBook_Size/wCOI overlays; wholePop adds PC outlier CSVs)
  §3  FST book-size-within-gait (Trotter/Pacer)
  §4  filtered.LD_prune.het_stats.${rg}.het (via --read-freq)
  §5  PCA COI overlay using per-group .het
  §6  bcftools roh on group-subset VCF + L1/L2/L3 filter chain
  §7  per-base consensus ROH (wholePop alone, or gait + 3 book-size subs)
  §8  roh_summary_by_RG_L3_Froh.${rg}.txt (F_ROH summary)
  §9  per-group GRM + ROHRM + analysis_comparison at every $ROH_CUTOFFS cutoff
      → Inbreeding_Comparison.${rg}.csv + Pairwise_Differences.${rg}.csv
  §10 related_${rg} (gait-filtered first-degree pair list)
  §11 PCA pairwise Euclidean distance
  §12 Euclidean + KING merge
  §13 Euclidean vs KING kinship correlation plot
  §14 cross-method correlation plot (ROHRM vs KING vs PCA)
  §15 F_SNP / F_ROH / D_STD / D_ROH / ROH_sh correlation plots (--mode froh,
      --mode froh-cons, and four doubleAnn plots)
  ↓
genDiversity_aggregate.sh  (runs once)
  → twoGait / threeBooksize per-sample ROH_sh concatenations
  → F_ROH histograms + roh_high.csv + gait/book-size stratified F_ROH summaries
  → Froh-vs-ROHsh plots (wholePop / twoGait / threeBooksize)
  → merged-kin top-pair cross-reference
  ↓
Google Drive (rclone upload — interleaved, not a dedicated phase)
```

### Repository layout expectations

The pipeline creates and uses one top-level directory per run:

- `results_<timestamp>/` — per-run output directory. Contains the downloaded inputs (`SNPdata_iScan_Standardbred/`, `Miscellaneous_documents_standardbred/`) alongside pipeline outputs (`preprocess/`, `dedup/`, `inspect/`, `filtered/`, `LD_pruned/`, `divStats/`, `rep_ROHRM/`) and the run log (`run.log`). `<timestamp>` is `YYYYMMDD_HHMMSS` captured when the script starts; override by exporting `OUTPUT_DIR=<existing_dir>` before invocation to reuse or resume into a prior folder. Subdir names after the timestamped prefix are unchanged from the pre-refactor layout. Per-group work adds `rep_ROHRM/perGroup_${rg}/` working dirs and `preprocess/samples.${rg}.txt` / `sample_groups.tsv` artifacts.

`results_*/` is in `.gitignore`.

Each run writes its full stdout+stderr to `${OUTPUT_DIR}/run.log` via a `tee` + `exec` redirection set up in `genDiversity_common.sh`. The redirection is guarded by the `GENDIV_LOG_SETUP` env var so subscripts invoked by the wrapper inherit its pipe instead of piling on their own tee and double-writing every line.

The pipeline also expects two sibling repositories and a configured rclone remote to exist:

- `../Equine80select_remapper/results/matchingSNPs_binary_consistantMapping.equCab3_map` — EquCab3 remap table
- `../Horse_parentage_SNPs/equCab3/download/equCab3.fa` + `equCab3_genome.fa.fai` — reference genome
- `remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs` — rclone remote for input download and output upload

A fresh clone of `genDiv` alone will fail early with path errors.

### GPA per-group reference files

`per_group.sh` produces deliverables consumed by the downstream GPA report pipeline (`../GPA/create_popFiles.sh`). The authoritative contract — filenames, schemas, column indices, and the per-file table of "numerically equivalent to today vs genuinely new content" — lives in **`MIGRATION.md`**. Any change to a producer or a consumer must update MIGRATION.md in the same PR.

Quick summary of the eight per-group files:

| File | Location | Purpose |
|---|---|---|
| `pruned.${rg}.afreq` | `LD_pruned/` | group-specific allele frequencies (feeds `--read-freq`) |
| `filtered.LD_prune.het_stats.${rg}.het` | `divStats/` | F_SNP het reference, computed with group AF |
| `roh_summary_by_RG_L3_Froh.${rg}.txt` | `divStats/` | F_ROH reference (bcftools roh on group-subset VCF) |
| `roh.L3.consensus_25pct.merged.${rg}.smoothed.bed` | `divStats/` | Consensus ROH regions from the group's own ROH calls |
| `roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt` | `divStats/` | ROH_sh per sample against the group's consensus |
| `Inbreeding_Comparison.${rg}.csv` | `rep_ROHRM/roh_1Mb.Threshold_3SD/` | D_SNP (col idx 1) + D_ROH (col idx 2); header `IID,D_STD,D_ROH,Phenotype` |
| `Pairwise_Differences.${rg}.csv` | `rep_ROHRM/roh_1Mb.Threshold_3SD/` | G_SNP (col idx 4) + G_ROH (col idx 5); header `ID1,ID2,Pheno1,Pheno2,Kinship_Std,Kinship_ROH,Difference` |
| `sample_groups.tsv` | `preprocess/` | Global sample → primary-group map (one row per sample) |

Every file follows the strict `<stem>.${rg}.<ext>` convention — wholePop outputs carry `.wholePop.` like Trotter and Pacer carry their own tags. `scripts/validate_popRefs.sh --mode upstream --root $OUTPUT_DIR` enforces the filename / header / row-count contract at the end of each run.

## Key Parameters & Thresholds

All tunable values are centralized in the CONFIG block at the top of `genDiversity.sh` (search for `# CONFIG`). Highlights:

- **Genome & data**
  - Reference: EquCab3
  - SNP array: Equine80select (~71,548 post-QC SNPs)
  - Autosomal genome length used to normalize F_ROH
- **Reproducibility & parallelism**
  - `BEAGLE_SEED=12345` (phasing seed)
  - `nthreads=10`
- **Quality control**
  - `GENO_MISS=0.05`, `MAF=0.01`, `HWE_PVAL=1e-6`
  - `KING_CUTOFF=0.20` (first-degree relative cutoff)
- **Sex calling** (X-chromosome F-stat)
  - `SEX_MAX_FEMALE_XF=0.2`, `SEX_MIN_MALE_XF=0.8`
  - `PAR_END_BP=2063653` (X-chr pseudoautosomal boundary on EquCab3)
- **LD pruning**
  - `LD_WINDOW_KB=100`, `LD_R2=0.8`
- **ROH & consensus**
  - `PRIMARY_ROH_MB=1.0` (primary window used downstream)
  - `ROH_CUTOFFS="1.0 5.0 10.0"` (additional cutoffs evaluated)
  - `pct=25` (min % samples in ROH to define consensus region)
  - `CONSENSUS_MIN_MB=0.5`
  - `ROH_THRESHOLD_SD=3.0` (window SNP-count filter)
- **Directory layout**
  - `OUTPUT_DIR="results_$(date +%Y%m%d_%H%M%S)"` — fresh timestamped folder every run; override by exporting `OUTPUT_DIR=<existing_dir>` before invocation to reuse or resume into a prior folder.
- **Phasing:** BEAGLE (required before `ROHRM_Creator.py`)

## Running Individual Python Scripts

Python scripts are called from within `genDiversity.sh` but can be run standalone:

```bash
conda activate genDiv
python scripts/ROHRM_Creator.py <args>
python scripts/analysis_comparison.py <args>
```

Check the corresponding section in `genDiversity.sh` for the exact arguments passed to each script.
