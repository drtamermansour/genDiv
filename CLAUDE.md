# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a genomic diversity assessment pipeline for Standardbred horses (~500+ individuals). It analyzes SNP genotyping data to compute genetic diversity metrics (ROH, heterozygosity, relatedness, PCA) and generate population-level reference ranges for breeding programs.

See also:
- `README.md` — human-oriented setup recipe (conda/mamba environment).
- `REFACTORING.md` — open suggestions and TODOs collected while reviewing the codebase.

## Running the Pipeline

The pipeline is a thin wrapper that composes three subscripts:

```bash
conda activate genDiv
bash ./genDiversity.sh
```

The wrapper sources shared CONFIG/helpers, invokes shared preprocessing once, loops the per-group stage three times (wholePop/Trotter/Pacer), then runs cross-group aggregation. Each subscript is independently runnable, which is handy for debugging or refreshing one stage:

```bash
bash ./genDiversity_shared.sh                 # preprocessing (Sections 1–3 + whole-pop ROH)
bash ./genDiversity_per_group.sh wholePop     # per-group metrics (repeat for Trotter, Pacer)
bash ./genDiversity_per_group.sh Trotter
bash ./genDiversity_per_group.sh Pacer
```

## Environment Setup

Uses conda/mamba with a named environment `genDiv`. Full setup command lives in `README.md`; below is a summary focused on what the pipeline actually invokes at runtime.

**Runtime dependencies called by `genDiversity.sh`:**
- CLI tools: `plink`, `plink2`, `bcftools`, `beagle`, `rclone`
- Interpreters: `python` (with `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`), `Rscript` (with `ggplot2`, `gridExtra`, `viridis`, `reshape2`, `GGally`, `effsize`)

**In the env recipe but not currently invoked by the pipeline:** `gcta`, `bedtools`, `snakemake`. See `REFACTORING.md` §5 — decide whether to drop or wire in.

## Architecture

### Entry Point

The pipeline is split across four bash files at the repo root:

| File | Role | Invocation |
|---|---|---|
| `genDiversity.sh` | Thin wrapper: sources common, runs shared → loops per_group × 3 → runs aggregate | `bash genDiversity.sh` |
| `genDiversity_common.sh` | CONFIG, helpers (`log`, `run_python`, `run_r`, `upload`), `ERR` trap, log redirection (guarded by `GENDIV_LOG_SETUP` so subscripts don't double-log) | sourced |
| `genDiversity_shared.sh` | Whole-pop preprocessing that runs once: Sections 1–3 + whole-pop ROH calling + per-base consensus + effective-genome-length + whole-pop F_SNP het + `samples.{wholePop,Trotter,Pacer}.txt` | `bash genDiversity_shared.sh` |
| `genDiversity_per_group.sh` | Per-group metric stage; takes `$rg ∈ {wholePop, Trotter, Pacer}`, runs 3× | `bash genDiversity_per_group.sh <rg>` |

All scripts respect `OUTPUT_DIR` as an env override, so re-running into an existing folder (e.g., to refresh just one group) works the same way as a fresh timestamped run.

`set -eo pipefail` applies everywhere; the `ERR` trap reports the failing line number from whichever subscript errored.

Conceptual workflow phases (the split is organizational — execution order is the same as before):

1. **Data Download & Preprocessing** (in `shared.sh`) — rclone from Google Drive, PLINK ID updates, EquCab3 remapping, SNP deduplication, PLINK1→PLINK2 conversion.
2. **Data Exploration** (in `shared.sh`) — sex validation (X chr F-stats), PAR removal, HWE analysis.
3. **Final Filtering** (in `shared.sh`) — apply missingness/MAF/HWE thresholds to produce the clean dataset.
4. **ROH Analysis** (in `shared.sh`) — bcftools roh on the whole-pop phased VCF, L1/L2/L3 filter chain, per-base consensus regions across wholePop + gait + book-size subgroups, effective autosomal genome length. The book-size / twoGait / threeBooksize concatenations used by downstream plots are also produced here.
5. **Whole-pop F_SNP het** (in `shared.sh`) — `plink2 --het` on all samples, producing the whole-pop het file consumed by the per-group COI overlays.
6. **Per-group reference files** (in `per_group.sh`, run once per `$rg`) — PCA + overlays (wSex / wGait / wBook_Size / wCOI, with wSex and wGait wholePop-only), FST book-size-within-gait (Trotter/Pacer only), plus the GPA-proposal deliverables: `pruned.${rg}.afreq`, `filtered.LD_prune.het_stats.${rg}.het`, `roh_summary_by_RG_L3_Froh.${rg}.txt`, `Inbreeding_Comparison.${rg}.csv`, `Pairwise_Differences.${rg}.csv`.
7. **Whole-pop relatedness** (still inline in `genDiversity.sh` — Section 6 work, to be extracted later) — KING-robust kinship, IBS, PCA-based Euclidean distance, cross-method correlations.
8. **Upload** — `rclone` is invoked throughout each subscript; there's no single upload phase.

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
Google Drive (rclone download)
  ↓
genDiversity_shared.sh  (runs once)
  PLINK preprocessing → QC & filtering → LD pruning
  → bcftools roh on whole-pop phased VCF → L1/L2/L3 filter → per-base consensus
  → effective_autosomal_genome_length.txt
  → whole-pop plink2 --het (F_SNP seed)
  → preprocess/samples.{wholePop,Trotter,Pacer}.txt + sample_groups.tsv
  ↓
genDiversity_per_group.sh  (runs 3× — wholePop, Trotter, Pacer)
  plink2 --pca (+ wSex/wGait/wBook_Size/wCOI overlays, wholePop extras)
  → FST book-size-within-gait (Trotter/Pacer)
  → pruned.${rg}.afreq  → filtered.LD_prune.het_stats.${rg}.het (via --read-freq)
  → bcftools roh on group-subset VCF → roh_summary_by_RG_L3_Froh.${rg}.txt
  → plink2 --make-rel (group GRM) + ROHRM_Creator.py (group ROHRM)
  → analysis_comparison.py → Inbreeding_Comparison.${rg}.csv + Pairwise_Differences.${rg}.csv
  ↓
genDiversity.sh  (still inline — wholePop KING / IBS / Euclidean / correlation plots)
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

A fresh clone of `genDiv` alone will fail early with path errors. `REFACTORING.md` §7 tracks a suggestion to validate these at script start.

### GPA per-group reference files

The `per_group.sh` stage produces the deliverables consumed by the downstream GPA report pipeline (`../GPA/create_popFiles.sh`):

| File | Location | Purpose |
|---|---|---|
| `pruned.${rg}.afreq` | `LD_pruned/` | group-specific allele frequencies (feeds `--read-freq`) |
| `filtered.LD_prune.het_stats.${rg}.het` | `divStats/` | F_SNP het reference, computed with group AF |
| `roh_summary_by_RG_L3_Froh.${rg}.txt` | `divStats/` | F_ROH reference (bcftools roh on group-subset VCF) |
| `Inbreeding_Comparison.${rg}.csv` | `rep_ROHRM/roh_1Mb.Threshold_3SD/` | D_SNP (col idx 1) + D_ROH (col idx 2), IID,D_STD,D_ROH,Phenotype |
| `Pairwise_Differences.${rg}.csv` | `rep_ROHRM/roh_1Mb.Threshold_3SD/` | G_SNP (col idx 4) + G_ROH (col idx 5), ID1,ID2,Pheno1,Pheno2,Kinship_Std,Kinship_ROH,Difference |

GPA reads by fixed column index — any schema drift breaks the report silently. Keep `analysis_comparison.py` column orderings stable.

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
