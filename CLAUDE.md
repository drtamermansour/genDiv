# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a genomic diversity assessment pipeline for Standardbred horses (~500+ individuals). It analyzes SNP genotyping data to compute genetic diversity metrics (ROH, heterozygosity, relatedness, PCA) and generate population-level reference ranges for breeding programs.

See also:
- `README.md` — human-oriented setup and orientation: what a run produces, disk/runtime cost, how to run one stage.
- `environment.yml` — single source of truth for the conda environment; `scripts/check_env.py` proves it still covers the code.
- `MIGRATION.md` — the filename / schema contract between this pipeline and the downstream GPA report pipeline (`../GPA/`). Update it whenever a producer here or a consumer there changes.
- `scripts/validate_popRefs.sh` — dual-mode validator for the per-group reference files; runs at the end of this pipeline and at the end of GPA's `create_popFiles.sh` (both sides are wired up). The GPA copy has diverged by adding a `POPFILES_INVARIANT_SPECS` array — see `MIGRATION.md` §5 before porting changes across.
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

**`environment.yml` is the single source of truth.** Do not maintain a second dependency list here or in `README.md` — point at that file instead. It declares only direct dependencies (something a script invokes or imports), each annotated with its consumer, and leaves transitives to the solver.

```bash
mamba env create -f environment.yml && conda activate genDiv
```

Two invariants worth knowing before you touch dependencies:

- **`rclone` is deliberately not declared.** It comes from the cluster module system (`module load rclone` at `genDiversity_shared.sh:12` and `explore.sh:25`), not conda. `scripts/check_env.py` knows this via its `EXTERNALLY_PROVIDED` set.
- **Never regenerate `environment.yml` with `conda env export`.** It would record one machine's transitive closure and destroy the direct/transitive distinction the file exists to preserve.

**After adding or changing any script, run the drift check:**

```bash
python3 scripts/check_env.py            # exits 1 on any undeclared dependency
python3 scripts/check_env.py --strict   # also fails on declared-but-unused
```

This is not optional politeness — the failure it catches is silent. A package can be missing from the spec while the pipeline still runs, because something else pulled it in transitively; nothing breaks until the environment is built from scratch. That is exactly how `r-dplyr`, `r-tidyr`, and `tqdm` stayed undeclared for months.

If you add a dependency whose import name differs from its conda package name (e.g. `allel` → `scikit-allel`), add the mapping to `PACKAGE_ALIASES` in `scripts/check_env.py`, or the check will report a false positive.

`gcta` and `snakemake` are commented out in `environment.yml` — installed historically, invoked by nothing today.

## Architecture

### Entry Point

The pipeline is split across five bash files at the repo root:

| File | Role | Invocation |
|---|---|---|
| `genDiversity.sh` | Thin wrapper: sources common, runs shared → loops per_group × 3 → runs aggregate | `bash genDiversity.sh` |
| `genDiversity_common.sh` | CONFIG, helpers (`log`, `run_python`, `run_r`, `upload`), `ERR` trap, log redirection (guarded by `GENDIV_LOG_SETUP` so subscripts don't double-log) | sourced |
| `genDiversity_shared.sh` | Whole-pop preprocessing (Sections 1–3) + effective-genome-length + autosomes.genome + `samples.{wholePop,Trotter,Pacer}.txt` + `sample_groups.tsv` | `bash genDiversity_shared.sh` |
| `genDiversity_per_group.sh` | Per-group metric stage; takes `$rg ∈ {wholePop, Trotter, Pacer}`, runs 3× | `bash genDiversity_per_group.sh <rg>` |
| `genDiversity_aggregate.sh` | Cross-group, five jobs: `fst_stats.R` over all five FST summaries; twoGait / threeBooksize ROH_sh concatenations + Froh-vs-ROHsh plots; ROH_common cross-group figures + pairwise stats; ROH islands + gene annotation; the LD-decay *N_e* tool loop | `bash genDiversity_aggregate.sh` |

All scripts respect `OUTPUT_DIR` as an env override, so re-running into an existing folder (e.g., to refresh just one group) works the same way as a fresh timestamped run.

`set -eo pipefail` applies everywhere; the `ERR` trap reports the failing line number from whichever subscript errored.

Conceptual workflow phases:

1. **Data Download & Preprocessing** (`shared.sh`) — rclone from Google Drive, PLINK ID updates, EquCab3 remapping, SNP deduplication, PLINK1→PLINK2 conversion.
2. **Data Exploration** (`shared.sh`) — sex validation (X chr F-stats), PAR removal, HWE analysis.
3. **Final Filtering** (`shared.sh`) — apply missingness/MAF/HWE thresholds to produce the clean dataset; LD pruning.
4. **Whole-pop preprocessing tail** (`shared.sh`) — derives `effective_autosomal_genome_length.txt`, `autosomes.genome`, and per-group sample lists (`samples.${rg}.txt` + `sample_groups.tsv`) so `per_group.sh` can run.
5. **Per-group reference files** (`per_group.sh`, 3×) — sections mirroring the original pipeline order. Highlights: per-group afreq, PCA + overlays (wSex / wGait wholePop-only, wBook_Size / wCOI all groups), A_e and whole-pop FST (wholePop-only), FST book-size-within-gait, per-group F_SNP `.het`, per-group bcftools roh + L1/L2/L3 + consensus ROH (nested over book-size for Trotter/Pacer), per-group tabix `freqs.${rg}.tab.gz` emitted alongside the ROH call for downstream GPA `bcftools roh --AF-file` consumption, F_ROH summary, per-group GRM + ROHRM + analysis_comparison across all `$ROH_CUTOFFS`, per-group KING + IBS + `related.${rg}`, PCA pairwise Euclidean, cross-method correlation plot, and COI-vs-F_ROH / F_ROH-vs-D_ROH / F_SNP-vs-D_ROH / F_SNP-vs-F_ROH doubleAnn plots.
6. **Cross-group aggregation** (`aggregate.sh`) — `fst_stats.R` over the five FST summaries, twoGait and threeBooksize per-sample ROH_sh concatenations, Froh-vs-ROHsh plots iterating wholePop / twoGait / threeBooksize, ROH_common cross-group figures + pairwise subgroup stats, and ROH islands with gene annotation.
7. **Effective population size** (`aggregate.sh`) — the opt-out `scripts/ne/` loop runs each installed LD-decay *N_e* tool (GONE2, currentNe2, NeEstimator, SNeP) once per group into `divStats/ne/`. See "LD-decay effective population size" below; this is a full analysis stage, not a sub-step of phase 6.
8. **Upload** — `rclone` is invoked throughout each subscript; there's no single upload phase.

### Python Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `ROHRM_Creator.py` | Builds the ROH-based Relationship Matrix (ROHRM) — implements Howard et al. C++ logic in Python using window-based haplotype matching on phased VCF data |
| `analysis_comparison.py` | Compares ROHRM vs standard GRM; `RobustMatrixComparator` class |
| `summary_roh.py` | ROH summary statistics grouped by subpopulation. `-n/--n-numeric-cols N` (default 3) controls how many numeric columns are summarised; the column right after the numeric block is the grouping factor |
| `summary_nseg_bins.py` | Per-subgroup summary of L3 ROH segment counts binned by length: `NSEG_1to3`, `NSEG_3to5`, `NSEG_5to10`, `NSEG_more10` (Mb, half-open [low, high)). Takes the per-segment `roh.L3.${rg}.txt` and a factor TSV (IID in col 1, factor in last col, same shape as `roh.L3_Froh_gait.txt`) |
| `summary_grm_kinship.py` | Per-subgroup pairwise GRM kinship summary. Reads the wholePop PLINK2 `.rel` (square) + `.rel.id` and the `USTA_Diversity_Study.gait_bookSize` factor; writes a 9-row CSV of within-subgroup off-diagonal pairwise GRM values (wholePop + 2 gaits + 6 gait × book-size strata) plus a violin plot. Each row carries mean ± SD and `Pct_pairs_gt_0.20` — the percent of pairs above `--tail-threshold` (default 0.20, the KING first-degree QC cutoff). The tail frequency, not the per-stratum maximum, is the reportable upper-tail statistic: the maximum scales with the pair count, so the larger MEDIUM strata top out higher than HIGH despite HIGH carrying far denser tails. Uses the wholePop GRM (single common AF basis) so cross-gait subgroup comparison is on a single yardstick; the same VanRaden additive-genetic similarities feed PLINK2 `--pca`. |
| `pc_outlier_kinship.py` | Confirms that PC2/PC3/PC4 "hidden familial structure" corresponds to genuine high-kinship clusters. For each PC, identifies the top-N (default 10) individuals at each tail (most positive and most negative loadings) and compares the within-cluster mean off-diagonal pairwise GRM kinship to the cohort mean (≈ 0 by GRM construction). Emits `divStats/PC_outlier_kinship.csv` with one row per (PC, tail). |
| `roh_islands_annotate.py` | Identifies ROH islands per group from the per-window ROH frequency landscapes produced by `roh_common_landscape.py`, thresholds at the top-1% f_w (with absolute f_w ≥ 0.5 flagged as high-confidence), merges contiguous high-f_w windows allowing a small gap (≤ 2 windows by default), filters islands narrower than 500 kb, and annotates each island with overlapping Ensembl EquCab3 protein-coding genes. A curated horse-selection candidate-gene list (`scripts/horse_selection_candidates.tsv`) flags islands overlapping known selection loci (DMRT3, MSTN, LCORL/NCAPG, MC1R, KIT, ASIP, STX17, MITF, MEF2C, …). A cross-group consolidation pass `bedtools-merges` per-group island sets into unique regions and records which group(s) contributed. Emits `divStats/roh_islands/roh_islands.{wholePop,Trotter,Pacer}.csv` plus `divStats/roh_islands/roh_islands.consolidated.csv`. |
| `summary_het.py` | Heterozygosity summary stats (observed/expected, F-coefficients) |
| `roh_plot.py` | Scatter plots correlating F_ROH vs consensus ROH sharing |
| `roh_histograms.py` | Unified histogram script: use `--metric ratio` (ROH_shared/F_ROH) or `--metric shared` (Percent_of_Consensus_ROH) |
| `roh_common_landscape.py` | Builds the per-window population ROH frequency landscape (f_w = n_w/N) for one group — genome-wide + 4 length-class landscapes (`landscape.${rg}.{,1to3,3to5,5to10,more10}.tsv`). bedtools subprocess for all interval arithmetic. |
| `roh_common_individual.py` | Per-individual ROH_common scoring with **exact** leave-one-out: `mean over w in W_i of (n_w − 1)/(N − 1)`. Emits `roh_common.${rg}.tsv` with genome-wide + 4 class scores + window counts; NA when a length class is empty for an individual. |
| `roh_common_plot.py` | Multi-mode plotter for Figures 2–4 of the ROH_common manuscript: `manhattan` (raw line plot, no smoothing), `scatter`, `lengthbox`. |
| `roh_common_subgroup_summary.py` | 9-row mean +/- SD table across the 5 ROH_common metrics: wholePop (cohort-wide scoring) + Pacer/Trotter + 6 gait × book-size subgroups (within-gait scoring). Emits `divStats/roh_common/roh_common_subgroup_summary.csv`. |
| `roh_common_pairwise_stats.py` | Pairwise Mann-Whitney U + Cohen's d for ROH_common across the six gait × book-size subgroups, over 5 metrics (genome-wide + 4 length classes) — 15 unordered pairs × 5 metrics = 75 rows. `p_adj_bonferroni` is the Bonferroni FWER-adjusted p-value computed within each metric family (15 tests). Emits `divStats/roh_common/roh_common_pairwise_stats.tsv`. |
| `utils.py` | Shared constants (`BOOK_SIZE_ORDER`, `BOOK_SIZE_COLORS`) and `format_stats()` used by summary scripts |

### LD-decay effective population size (`scripts/ne/`)

Each LD-decay-based *N_e* tool lives in its own wrapper under `scripts/ne/`. The convention is **modular and opt-out**: `genDiversity_aggregate.sh` iterates `gone2`, `currentne2`, `neestimator`, `snep` and invokes any wrapper that exists and is executable. To remove a tool: delete (or `chmod -x`) its `run_*.sh` and `install_*.sh` files. No edits to `genDiversity_aggregate.sh` required.

GONE2 gets the **unpruned** post-QC set (its estimator needs the full LD spectrum); the other three get the **LD-pruned** set, matching McGivney 2020 / Manunza 2025 practice. Every wrapper runs once per group (wholePop / Trotter / Pacer) and emits a standardised summary CSV with columns `Group, Method, Generations_ago, Ne, CI_level, CI_low, CI_high`. `Generations_ago` is an integer generation for the trajectory tools and the literal `contemporary` for the single-point tools.

| Wrapper | Tool | Reference | Input | Output |
|---|---|---|---|---|
| `scripts/ne/run_gone2.sh` (+ `install_gone2.sh`) | **GONE2 v1.0.2** (Santiago, Köpke & Caballero 2025; doi:10.1038/s41467-025-61378-w) | per-generation *N_e* trajectory from a single SNP sample using the full LD spectrum | post-QC PLINK 1 `.bed/.bim/.fam` set, **not** LD-pruned (full LD spectrum required); MAF 0.05; constant rec rate 1.16 cM/Mb (Beeson et al. 2020) | `divStats/ne/gone2/<group>/<group>_GONE2_{Ne,STATS,d2}` + `divStats/ne/gone2/Ne_gone2_summary.csv` |
| `scripts/ne/run_currentne2.sh` (+ `install_currentne2.sh`) | **currentNe2** (Santiago, Köpke & Caballero 2025) | contemporary single-point *N_e* from LD between mostly-unlinked loci; the faster analogue to NeEstimator recommended by the GONE2 paper | LD-pruned PLINK set; rec rate 1.16 cM/Mb | `divStats/ne/currentne2/<group>/input.<group>_currentNe2_OUTPUT.txt` + `Ne_currentne2_summary.csv`. `CI_level=90` — currentNe2 emits 50% / 90% CIs natively, not 95% |
| `scripts/ne/run_neestimator.sh` (+ `install_neestimator.sh`) | **NeEstimator v2.x** (Do et al. 2014) | contemporary point *N_e* via the LD method; method-matched to McGivney 2020's Thoroughbred *N_e* = 330 | LD-pruned PLINK set converted to GENEPOP; P_Crit 0.02, random mating, Waples (2006) bias correction, jackknife CIs (Manunza 2025 livestock recommendations) | `divStats/ne/neestimator/<group>/input.<group>Ne.txt` + `Ne_neestimator_summary.csv`. `CI_level=95_jackknife` when the jackknife block parses |
| `scripts/ne/run_snep.sh` (+ `install_snep.sh`) | **SNeP v1.11** (Barbato et al. 2015; doi:10.3389/fgene.2015.00109) | per-generation *N_e* trajectory from LD decay via the Sved–Feldman (1971) approximation; the dominant tool in 2023–2025 livestock literature | LD-pruned PLINK set; MAF 0.05, rec rate 1.16e-8 M/bp, `-itemsTH 500` minimum SNP pairs per distance bin | `divStats/ne/snep/<group>/*.{NeAll,LDAll}` + `Ne_snep_summary.csv` |
| `scripts/ne/plot_ne_trajectory.py` | — | overlays one or more summary CSVs on a shared *N_e*-vs-generations axis: log-y, a split x-axis at gen 15, a shaded gen 1–4 GONE artifact band (Novo et al. 2023), Standardbred event lines (2009 USTA cap, 1973 studbook closure, ~1872 breed founding), and a calendar-year secondary axis at G = 11.0 yr/gen (Waples 2026) | any standardised `Ne_<tool>_summary.csv` | `<out>.png` + `<out>.pdf` |

**GONE2 metapopulation sensitivity run.** `run_gone2.sh --metapopulation` passes GONE2's `-x` flag (sample drawn from a metapopulation of equal-sized subpopulations) instead of the default panmixia assumption. Pair it with `--output-subdir gone2_x --method-label GONE2_x` so the sensitivity outputs land in `divStats/ne/gone2_x/` and stay distinguishable in the summary CSV's `Method` column rather than clobbering the panmixia run. Invoked manually, not by the `aggregate.sh` loop.

**Interpretation caveat.** GONE2's gen 1–4 values are a single repeated number — a block-estimator saturation artifact, not signal. Do not read that plateau as a biological cross-group difference; start at gen 5.

All Python scripts use `argparse`; run with `--help` to see usage.

### R Scripts (`scripts/`)

Grouped by role. All take positional arguments via `commandArgs(trailingOnly = TRUE)`; grep the owning stage (`genDiversity_per_group.sh` / `genDiversity_aggregate.sh`) for a real call site.

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

**Statistics / plotting scripts** (take the run directory as their first positional argument, then resolve their inputs beneath it):

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
  → effective_autosomal_genome_length.txt + autosomes.genome
  → preprocess/samples.{wholePop,Trotter,Pacer}.txt + sample_groups.tsv
  ↓
genDiversity_per_group.sh  (runs 3× — wholePop, Trotter, Pacer)
  §1  pruned.${rg}.afreq
  §2  plink2 --pca (+ wSex/wGait/wBook_Size/wCOI overlays; wholePop adds PC outlier CSVs)
  §3  FST book-size-within-gait (Trotter/Pacer)
  §4  filtered.LD_prune.het_stats.${rg}.het (via --read-freq)
  §5  PCA COI overlay using per-group .het
  §6  bcftools roh on group-subset VCF + L1/L2/L3 filter chain
      (also emits freqs.${rg}.tab.gz + .tbi for downstream GPA AF-file consumption)
  §7  per-base consensus ROH (wholePop alone, or gait + 3 book-size subs)
  §7.5 ROH_common landscape (genome-wide + 4 length-class f_w maps for the
       main group) + per-individual ROH_common with exact LOO;
       writes divStats/roh_common/
  §8  roh_summary_by_RG_L3_Froh.${rg}.txt (F_ROH summary)
  §8b wholePop-only: gait / gait_bookSize F_ROH summary CSVs +
      roh.L3_NSEGbins_{gait,gait_bookSize}.sumStats.csv (NSEG length-bin summaries)
  §9  per-group GRM + ROHRM + analysis_comparison at every $ROH_CUTOFFS cutoff
      → Inbreeding_Comparison.${rg}.csv + Pairwise_Differences.${rg}.csv
      §9c (wholePop only) per-subgroup pairwise GRM kinship summary →
      divStats/GRM_kinship_by_subgroup.{csv,violin.png}
      §9d (wholePop only) PC2-PC4 outlier-cluster kinship confirmation →
      divStats/PC_outlier_kinship.csv
  §10 per-group KING + IBS + related.${rg} (per-group --make-king-table)
  §11 PCA pairwise Euclidean distance
      (§12 and §13 do not exist — numbering is inherited from the
       pre-refactor script and was left alone to keep §14/§15 stable)
  §14 cross-method correlation plot (ROHRM vs KING vs PCA)
  §15 F_SNP / F_ROH / D_STD / D_ROH / ROH_sh correlation plots (--mode froh,
      --mode froh-cons, and four doubleAnn plots)
  ↓
genDiversity_aggregate.sh  (runs once)
  → fst_stats.R (reads all 5 FST summaries: 3 whole-pop + 2 per-gait book-size)
  → twoGait / threeBooksize per-sample ROH_sh concatenations
  → Froh-vs-ROHsh plots (wholePop / twoGait / threeBooksize)
  → ROH_common: twoGait / threeBooksize per-sample score concatenations +
                Figures 2-4 (Manhattan landscape, F_ROH vs ROH_common
                scatter, length-stratified boxplots)
  → ROH islands: per-group top-1% f_w islands + cross-group
                 consolidation; gene annotation against Ensembl EquCab3
                 GTF + curated horse-selection candidate-gene list →
                 divStats/roh_islands/
  → N_e: for tool in gone2 currentne2 neestimator snep, run
         scripts/ne/run_<tool>.sh if executable, once per group →
         divStats/ne/<tool>/Ne_<tool>_summary.csv (+ per-group raw
         outputs). Skipped silently for any tool whose wrapper is
         absent or non-executable.
  ↓
Google Drive (rclone upload — interleaved, not a dedicated phase)
```

### Repository layout expectations

The pipeline creates and uses one top-level directory per run:

- `results_<timestamp>/` — per-run output directory. Contains the downloaded inputs (`SNPdata_iScan_Standardbred/`, `Miscellaneous_documents_standardbred/`) alongside pipeline outputs (`preprocess/`, `dedup/`, `inspect/`, `filtered/`, `LD_pruned/`, `divStats/`, `rep_ROHRM/`) and the run log (`run.log`). `<timestamp>` is `YYYYMMDD_HHMMSS` captured when the script starts; override by exporting `OUTPUT_DIR=<existing_dir>` before invocation to reuse or resume into a prior folder. Subdir names after the timestamped prefix are unchanged from the pre-refactor layout. Per-group work adds `rep_ROHRM/perGroup_${rg}/` working dirs and `preprocess/samples.${rg}.txt` / `sample_groups.tsv` artifacts.

`results_*/` is in `.gitignore`. So are `manuscript/` and `genDiv_manuscript/` — local-only writing directories that are not part of the pipeline.

### One-off exploratory analyses (`explore.sh`, `explore/`)

`explore.sh` is an opt-in driver for exploratory questions that do not belong in the main pipeline. It sources `genDiversity_common.sh` for the `log` / `run_python` / `upload` helpers, reads an **existing** run directory, writes to `${OUTPUT_DIR}/explore/`, and rclone-uploads each deliverable to `$GDRIVE_BASE/explore/<task_name>/`.

```bash
OUTPUT_DIR=results_<timestamp> bash explore.sh              # all tasks
OUTPUT_DIR=results_<timestamp> bash explore.sh <task_name>  # one task
```

Each task is a self-contained block calling a script in `explore/`, with a header listing its question, inputs, local outputs, and GDrive destination. New tasks follow that convention and must upload their deliverables. Nothing here runs as part of `genDiversity.sh`.

| Script | Purpose |
|---|---|
| `explore/relationship_comparison_colored.py` | Re-draws the Standard-GRM vs ROH-GRM pair-kinship scatter from `Robust_Matrix_Comparison_Enhanced.wholePop.png`, recoloured by pair-level gait (3 categories) and by pair-level book size (6 categories), one panel per ROH cutoff (1 / 5 / 10 Mb) with axes shared within a figure. The book-size view is emitted three times: all six categories, plus HIGH-\* and LOW/MEDIUM-\* subsets that keep the full-view palette but re-fit their axes. |

Each run writes its full stdout+stderr to `${OUTPUT_DIR}/run.log` via a `tee` + `exec` redirection set up in `genDiversity_common.sh`. The redirection is guarded by the `GENDIV_LOG_SETUP` env var so subscripts invoked by the wrapper inherit its pipe instead of piling on their own tee and double-writing every line.

The pipeline also expects two sibling repositories and a configured rclone remote to exist:

- `../InfiniTier/results_E80selv2_to_equCab3noAlt_genDiv/qc/Equine80select_v2_1_HTS_20143333_B1_UCD_allele_map_equCab3noAlt.tsv` — EquCab3 remap table (`$equCab3_map`)
- `../Horse_parentage_SNPs/equCab3/download/equCab3.fa` (`$ref`) + `../Horse_parentage_SNPs/equCab3/equCab3_genome.fa.fai` (`$reference_fai`) — reference genome
- `remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs` (`$GDRIVE_BASE`) — rclone remote for input download and output upload

All four are set at the top of the CONFIG block in `genDiversity_common.sh`; they are the first things to edit when moving the pipeline to a new machine. A fresh clone of `genDiv` alone will fail early with path errors.

### GPA per-group reference files

`per_group.sh` produces deliverables consumed by the downstream GPA report pipeline (`../GPA/create_popFiles.sh`). The authoritative contract — filenames, schemas, column indices, and the per-file table of "numerically equivalent to today vs genuinely new content" — lives in **`MIGRATION.md`**. Any change to a producer or a consumer must update MIGRATION.md in the same PR.

Quick summary of the nine per-group files:

| File | Location | Purpose |
|---|---|---|
| `pruned.${rg}.afreq` | `LD_pruned/` | group-specific allele frequencies on the LD-pruned SNP set (feeds PLINK2 `--read-freq`) |
| `freqs.${rg}.tab.gz` (+ `.tbi`) | `divStats/` | tabix-indexed `CHROM POS REF,ALT AF` on the filtered SNP set; consumed by downstream GPA's `bcftools roh --AF-file` on 1–2 animal mate-pair VCFs |
| `filtered.LD_prune.het_stats.${rg}.het` | `divStats/` | F_SNP het reference, computed with group AF |
| `roh_summary_by_RG_L3_Froh.${rg}.txt` | `divStats/` | F_ROH reference (bcftools roh on group-subset VCF) |
| `roh.L3.consensus_25pct.merged.${rg}.smoothed.bed` | `divStats/` | Consensus ROH regions from the group's own ROH calls |
| `roh.L3.perSample_intersect_${rg}_consensus_25pct.summary.txt` | `divStats/` | ROH_sh per sample against the group's consensus |
| `Inbreeding_Comparison.${rg}.csv` | `rep_ROHRM/roh_1Mb.Threshold_3SD/` | D_SNP (col idx 1) + D_ROH (col idx 2); header `IID,D_STD,D_ROH,Phenotype` |
| `Pairwise_Differences.${rg}.csv` | `rep_ROHRM/roh_1Mb.Threshold_3SD/` | G_SNP (col idx 4) + G_ROH (col idx 5); header `ID1,ID2,Pheno1,Pheno2,Kinship_Std,Kinship_ROH,Difference` |
| `sample_groups.tsv` | `preprocess/` | Global sample → primary-group map (one row per sample) |

Every file follows the strict `<stem>.${rg}.<ext>` convention — wholePop outputs carry `.wholePop.` like Trotter and Pacer carry their own tags. `scripts/validate_popRefs.sh --mode upstream --root $OUTPUT_DIR` enforces the filename / header / row-count contract at the end of each run.

## Key Parameters & Thresholds

All tunable values are centralized in the CONFIG block at the top of **`genDiversity_common.sh`** (search for `# CONFIG`), which every subscript sources. Highlights:

- **External paths** (edit these first on a new machine)
  - `equCab3_map`, `ref`, `reference_fai` — sibling-repo inputs; see "Repository layout expectations"
  - `GDRIVE_BASE` — rclone remote for both download and upload
- **Genome & data**
  - Reference: EquCab3
  - SNP array: Equine80select v2.1 — 76,841 assayed → 71,548 after SNP deduplication (`genDiversity_shared.sh:177`) → ~58k after the missingness / MAF / HWE filters. Use the ~58k figure for anything downstream of Section 3; `MIGRATION.md` §3 records the exact per-run counts.
  - Autosomal genome length used to normalize F_ROH — derived at runtime into `effective_autosomal_genome_length.txt`, not a CONFIG constant
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
- **ROH_common (continuous population-autozygosity metric)**
  - `ROH_COMMON_WINDOW_KB=100` (fixed window size for the f_w landscape)
  - `ROH_COMMON_LENGTH_BINS="1,3,5,10"` (ROH length-class edges in Mb; aligned with `summary_nseg_bins.py` bins)
- **Directory layout**
  - `OUTPUT_DIR="results_$(date +%Y%m%d_%H%M%S)"` — fresh timestamped folder every run; override by exporting `OUTPUT_DIR=<existing_dir>` before invocation to reuse or resume into a prior folder.
- **Phasing:** BEAGLE (required before `ROHRM_Creator.py`)

## Running Individual Python Scripts

Every Python and R script can be run standalone against a finished `OUTPUT_DIR`:

```bash
conda activate genDiv
python scripts/ROHRM_Creator.py --help
python scripts/analysis_comparison.py --help
```

`genDiversity.sh` itself is a 26-line wrapper and calls no scripts directly. To find a real call site with its arguments, grep the stage that owns it:

```bash
grep -n "summary_grm_kinship" genDiversity_per_group.sh   # per-group stage
grep -n "roh_common_pairwise_stats" genDiversity_aggregate.sh   # cross-group stage
```

All Python scripts use `argparse`, so `--help` is authoritative for the interface; the call site is what tells you which files the pipeline actually feeds in.
