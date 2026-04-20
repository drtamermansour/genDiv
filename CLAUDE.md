# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a genomic diversity assessment pipeline for Standardbred horses (~500+ individuals). It analyzes SNP genotyping data to compute genetic diversity metrics (ROH, heterozygosity, relatedness, PCA) and generate population-level reference ranges for breeding programs.

See also:
- `README.md` — human-oriented setup recipe (conda/mamba environment).
- `REFACTORING.md` — open suggestions and TODOs collected while reviewing the codebase.

## Running the Pipeline

The entire pipeline is orchestrated by a single bash script:

```bash
conda activate genDiv
bash ./genDiversity.sh
```

The script is not designed to be run in sections interactively — individual sections are commented/uncommented as needed during development.

## Environment Setup

Uses conda/mamba with a named environment `genDiv`. Full setup command lives in `README.md`; below is a summary focused on what the pipeline actually invokes at runtime.

**Runtime dependencies called by `genDiversity.sh`:**
- CLI tools: `plink`, `plink2`, `bcftools`, `beagle`, `rclone`
- Interpreters: `python` (with `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`), `Rscript` (with `ggplot2`, `gridExtra`, `viridis`, `reshape2`, `GGally`, `effsize`)

**In the env recipe but not currently invoked by the pipeline:** `gcta`, `bedtools`, `snakemake`. See `REFACTORING.md` §5 — decide whether to drop or wire in.

## Architecture

### Entry Point

**`genDiversity.sh`** — monolithic bash pipeline with `set -eo pipefail`, an `ERR` trap that reports the failing line number, a CONFIG section at lines 4–35 (thresholds, paths, seeds centralized there), and helper functions `log()`, `run_python()`, `run_r()`, `upload()` at lines 40–50.

Structurally the script is **not** 8 cleanly numbered sections. It has 5 `log "Section …"` markers (`genDiversity.sh:60, 135, 366, 481, 1145`) plus many inline subsections. Conceptually the workflow phases are:

1. **Data Download & Preprocessing** — rclone from Google Drive, PLINK ID updates, EquCab3 remapping, SNP deduplication, PLINK1→PLINK2 conversion
2. **Data Exploration** — sex validation (X chr F-stats), PAR removal, heterozygosity, allele frequency, HWE analysis
3. **Final Filtering** — apply missingness/MAF/HWE thresholds to produce clean dataset
4. **ROH Analysis** — call ROH with PLINK2 (`--hom`, min 1Mb), compute consensus ROH regions (≥25% samples, min 500kb), F_ROH inbreeding coefficients, stratify by gait (Pacer/Trotter) and book size (HIGH/MEDIUM/LOW)
5. **Nucleotide Diversity** — pi calculations. **Under development** — this section is a stub/comments about SNP-array ascertainment bias. Pending decision (see `REFACTORING.md` §1): complete an array-aware π implementation, or remove.
6. **Genomic Relatedness** — standard GRM (vanRaden via PLINK2) and novel ROH-based Relatedness Matrix (ROHRM)
7. **Additional Relatedness** — KING-robust kinship, IBS, PCA-based Euclidean distance, cross-method correlations
8. **Upload** — results to Google Drive via rclone (interleaved throughout, not a dedicated final section)

### Python Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `ROHRM_Creator.py` | Builds the ROH-based Relationship Matrix (ROHRM) — implements Howard et al. C++ logic in Python using window-based haplotype matching on phased VCF data |
| `ROH_analysis.py` | ROH region identification and statistics from VCF files |
| `analysis_comparison.py` | Compares ROHRM vs standard GRM; `RobustMatrixComparator` class |
| `summary_roh.py` | ROH summary statistics grouped by subpopulation (3 numeric columns, factor at col 4). Used at `genDiversity.sh:950` |
| `summary_roh_v2.py` | Variant of `summary_roh.py` accepting a 4th numeric column and inlining `format_stats` instead of importing from `utils.py`. Used at `genDiversity.sh:1101, 1108`. Consolidation candidate — see `REFACTORING.md` §11 |
| `summary_het.py` | Heterozygosity summary stats (observed/expected, F-coefficients) |
| `roh_plot.py` | Scatter plots correlating F_ROH vs consensus ROH sharing |
| `roh_histograms.py` | Unified histogram script: use `--metric ratio` (ROH_shared/F_ROH) or `--metric shared` (Percent_of_Consensus_ROH) |
| `utils.py` | Shared constants (`BOOK_SIZE_ORDER`, `BOOK_SIZE_COLORS`) and `format_stats()` used by summary scripts |

Scripts in `scripts/sandbox/` are experimental variants not used in the main pipeline.

All Python scripts use `argparse`; run with `--help` to see usage.

### R Scripts (`scripts/`)

Grouped by role. Each script's `genDiversity.sh` call site in parentheses.

**PCA:**

| Script | Purpose |
|--------|---------|
| `pca_plots.R` | Unified PCA plots: positional args `<prefix> <eigenvec_suffix> <color_col> <output.png> [n_pcs=3] [color_type=factor\|numeric]` |

**Correlation heatmaps & scatterplot matrices** — four versions coexist, each wired in at a different pipeline stage with different input sets. Consolidation is on the backlog (see `REFACTORING.md` §2).

| Script | Used at | Purpose |
|--------|---------|---------|
| `correlation_plot_multiway.R` | `genDiversity.sh:796, 821, 846, 964` | Initial version: `<homo_file> <het_file> <out_prefix>` |
| `correlation_plot_multiway_v2.R` | `genDiversity.sh:1263` | Adds F_ROH stats input |
| `correlation_plot_multiway_v2e.R` | `genDiversity.sh:1270` | v2 + consensus-share input |
| `correlation_plot_multiway_v3.R` | `genDiversity.sh:1427` | Adds KING-IBS and Euclidean-distance inputs |

**Scatter plots:**

| Script | Used at | Purpose |
|--------|---------|---------|
| `plot_correlation.R` | `genDiversity.sh:1335, 1413` | Scatter of two columns with single-column color |
| `plot_correlation_withColorsAndShapes.R` | `genDiversity.sh:1281–1302` | Scatter with both color and shape mappings |
| `plot_correlation_withColors.R` | only in commented-out call at `genDiversity.sh:1276` | Currently unused — see `REFACTORING.md` §3 |

**Statistics / plotting scripts** (hard-coded input paths, no CLI args — pipeline invokes them with no arguments; outputs captured by stdout redirects where applicable):

| Script | Used at | Purpose |
|--------|---------|---------|
| `fst_stats.R` | `genDiversity.sh:690` | FST statistics across subpopulations; output captured to `divStats/fst_stats.txt` |
| `plot_Ae.R` | `genDiversity.sh:654` | Violin plots of effective allele number (Ae) by gait and book size |
| `effAllele_stats.R` | `genDiversity.sh:653` | Effect size analysis for Ae; output captured to `divStats/effAllele_stats.txt` |

### Data Flow

```
Google Drive (rclone)
  → PLINK preprocessing (bash)
  → PLINK2 QC & filtering
  → PLINK2 ROH calling → ROHRM_Creator.py (phased VCF required)
  → PLINK2 GRM → analysis_comparison.py
  → R scripts (visualization)
  → Google Drive (rclone upload)
```

### Repository layout expectations

`genDiversity.sh` expects two sibling repositories and a configured rclone remote to exist:

- `../Equine80select_remapper/results/matchingSNPs_binary_consistantMapping.equCab3_map` — EquCab3 remap table
- `../Horse_parentage_SNPs/equCab3/download/equCab3.fa` + `equCab3_genome.fa.fai` — reference genome
- `remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs` — rclone remote for input download and output upload

A fresh clone of `genDiv` alone will fail early with path errors. `REFACTORING.md` §7 tracks a suggestion to validate these at script start.

## Key Parameters & Thresholds

All tunable values are centralized in the CONFIG block at `genDiversity.sh:4-35`. Highlights:

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
- **Phasing:** BEAGLE (required before `ROHRM_Creator.py`)

## Running Individual Python Scripts

Python scripts are called from within `genDiversity.sh` but can be run standalone:

```bash
conda activate genDiv
python scripts/ROHRM_Creator.py <args>
python scripts/analysis_comparison.py <args>
```

Check the corresponding section in `genDiversity.sh` for the exact arguments passed to each script.
