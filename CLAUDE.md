# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a genomic diversity assessment pipeline for Standardbred horses (~500+ individuals). It analyzes SNP genotyping data to compute genetic diversity metrics (ROH, heterozygosity, relatedness, PCA) and generate population-level reference ranges for breeding programs.

## Running the Pipeline

The entire pipeline is orchestrated by a single bash script:

```bash
conda activate genDiv
bash ./genDiversity.sh
```

The script is not designed to be run in sections interactively — individual sections are commented/uncommented as needed during development.

## Environment Setup

Uses conda/mamba with a named environment `genDiv`:

```bash
mamba create -n genDiv -c conda-forge -c bioconda \
  r-base=4.5.2 plink plink2 bcftools gcta beagle bedtools rclone \
  python numpy pandas scipy matplotlib seaborn
```

## Architecture

### Entry Point

**`genDiversity.sh`** — monolithic bash pipeline with `set -eo pipefail` and a CONFIG section at the top (all thresholds, paths, and constants centralized there). Helper functions `log()`, `run_python()`, `run_r()`, and `upload()` are defined after the CONFIG block. 8 sequential sections:

1. **Data Download & Preprocessing** — rclone from Google Drive, PLINK ID updates, EquCab3 remapping, SNP deduplication, PLINK1→PLINK2 conversion
2. **Data Exploration** — sex validation (X chr F-stats), PAR removal, heterozygosity, allele frequency, HWE analysis
3. **Final Filtering** — apply missingness/MAF/HWE thresholds to produce clean dataset
4. **ROH Analysis** — call ROH with PLINK2 (`--hom`, min 1Mb), compute consensus ROH regions (≥25% samples, min 500kb), F_ROH inbreeding coefficients, stratify by gait (Pacer/Trotter) and book size (HIGH/MEDIUM/LOW)
5. **Nucleotide Diversity** — pi calculations (under development)
6. **Genomic Relatedness** — standard GRM (vanRaden via PLINK2) and novel ROH-based Relatedness Matrix (ROHRM)
7. **Additional Relatedness** — KING-robust kinship, IBS, PCA-based Euclidean distance, cross-method correlations
8. **Upload** — results to Google Drive via rclone

### Python Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `ROHRM_Creator.py` | Builds the ROH-based Relationship Matrix (ROHRM) — implements Howard et al. C++ logic in Python using window-based haplotype matching on phased VCF data |
| `ROH_analysis.py` | ROH region identification and statistics from VCF files |
| `analysis_comparison.py` | Compares ROHRM vs standard GRM; `RobustMatrixComparator` class |
| `summary_roh.py` | ROH summary statistics grouped by subpopulation |
| `summary_het.py` | Heterozygosity summary stats (observed/expected, F-coefficients) |
| `roh_plot.py` | Scatter plots correlating F_ROH vs consensus ROH sharing |
| `roh_histograms.py` | Unified histogram script: use `--metric ratio` (ROH_shared/F_ROH) or `--metric shared` (Percent_of_Consensus_ROH) |
| `utils.py` | Shared constants (`BOOK_SIZE_ORDER`, `BOOK_SIZE_COLORS`) and `format_stats()` used by summary scripts |

Scripts in `scripts/sandbox/` are experimental variants not used in the main pipeline.

All Python scripts use `argparse`; run with `--help` to see usage.

### R Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `pca_plots.R` | Unified PCA plots: positional args are `<prefix> <eigenvec_suffix> <color_col> <output.png> [n_pcs=3] [color_type=factor\|numeric]` |
| `correlation_plot_multiway*.R` | Pairwise correlation heatmaps and scatterplot matrices across inbreeding metrics |
| `fst_stats.R` | FST statistics across subpopulations |
| `plot_Ae.R` | Violin plots of effective allele number (Ae) by gait and book size |

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

## Key Parameters & Thresholds

- **Genome reference:** EquCab3
- **SNP array:** Equine80select (~71,548 post-QC SNPs)
- **ROH minimum size:** 1 Mb
- **Consensus ROH:** ≥25% samples in ROH, min 500 kb
- **Autosomal genome length:** used to normalize F_ROH
- **Phasing:** BEAGLE (required before ROHRM_Creator.py)

## Running Individual Python Scripts

Python scripts are called from within `genDiversity.sh` but can be run standalone:

```bash
conda activate genDiv
python scripts/ROHRM_Creator.py <args>
python scripts/analysis_comparison.py <args>
```

Check the corresponding section in `genDiversity.sh` for the exact arguments passed to each script.
