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
Google Drive (rclone)
  → PLINK preprocessing (bash)
  → PLINK2 QC & filtering
  → PLINK2 ROH calling → ROHRM_Creator.py (phased VCF required)
  → PLINK2 GRM → analysis_comparison.py
  → R scripts (visualization)
  → Google Drive (rclone upload)
```

### Repository layout expectations

`genDiversity.sh` creates and uses two top-level directories:

- `input_data/` — `SNPdata_iScan_Standardbred/` (downloaded genotypes) and `Miscellaneous_documents_standardbred/` (downloaded metadata).
- `results_<timestamp>/` — per-run output directory. `<timestamp>` is `YYYYMMDD_HHMMSS` captured when the script starts. Contains `preprocess/`, `dedup/`, `inspect/`, `filtered/`, `LD_pruned/`, `divStats/`, `rep_ROHRM/`, and the run log (`run.log`). Subdir names after the timestamped prefix are unchanged from the pre-refactor layout.

`input_data/` and `results_*/` are in `.gitignore`.

Each run writes its full stdout+stderr to `${OUTPUT_DIR}/run.log` via a `tee` + `exec` redirection set near the top of the script.

`genDiversity.sh` also expects two sibling repositories and a configured rclone remote to exist:

- `../Equine80select_remapper/results/matchingSNPs_binary_consistantMapping.equCab3_map` — EquCab3 remap table
- `../Horse_parentage_SNPs/equCab3/download/equCab3.fa` + `equCab3_genome.fa.fai` — reference genome
- `remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs` — rclone remote for input download and output upload

A fresh clone of `genDiv` alone will fail early with path errors. `REFACTORING.md` §7 tracks a suggestion to validate these at script start.

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
  - `INPUT_DIR="input_data"`
  - `OUTPUT_DIR="results_$(date +%Y%m%d_%H%M%S)"` — fresh timestamped folder every run
- **Phasing:** BEAGLE (required before `ROHRM_Creator.py`)

## Running Individual Python Scripts

Python scripts are called from within `genDiversity.sh` but can be run standalone:

```bash
conda activate genDiv
python scripts/ROHRM_Creator.py <args>
python scripts/analysis_comparison.py <args>
```

Check the corresponding section in `genDiversity.sh` for the exact arguments passed to each script.
