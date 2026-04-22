# genDiv

This pipeline makes use of 500+ horses' genotyping data to assess genetic diversity in the population and create population-level reference ranges for future breeding programs.

## Environment

Create a named conda/mamba environment (`genDiv`) containing every tool and package the pipeline invokes at runtime:

```bash
## (optional) update conda itself
conda update conda

## create the env
mamba create -n genDiv \
  conda-forge::r-base=4.5.2 \
  conda-forge::r-ggplot2=4.0.1 \
  conda-forge::r-gridextra=2.3 \
  conda-forge::r-qqman=0.1.9 \
  conda-forge::r-viridis=0.6.5 \
  conda-forge::r-reshape2=1.4.5 \
  conda-forge::r-ggally=2.4.0 \
  conda-forge::r-effsize=0.8.1 \
  conda-forge::openpyxl \
  conda-forge::scikit-allel \
  conda-forge::numpy conda-forge::pandas conda-forge::scipy \
  conda-forge::matplotlib conda-forge::seaborn \
  bioconda::plink bioconda::plink2 bioconda::bcftools bioconda::bedtools \
  bioconda::beagle bioconda::rclone

## Optional / not currently invoked by genDiversity.sh — keep if you plan to wire them in:
##   bioconda::gcta bioconda::snakemake=9.16.3

conda activate genDiv
```

## Running

The pipeline expects sibling repositories next to this one (`../Equine80select_remapper/` and `../Horse_parentage_SNPs/equCab3/`) and a configured rclone remote (`remote_UCDavis_GoogleDr`). See `CLAUDE.md` → "Repository layout expectations" for details.

```bash
git clone git@github.com:drtamermansour/genDiv.git
cd genDiv
conda activate genDiv
bash ./genDiversity.sh
```

The pipeline is split across five bash files at the repo root (`genDiversity.sh` is the wrapper; `genDiversity_common.sh` holds shared CONFIG/helpers; `genDiversity_shared.sh`, `genDiversity_per_group.sh`, and `genDiversity_aggregate.sh` do the work). Each subscript is independently runnable, useful for refreshing just one stage:

```bash
OUTPUT_DIR="Path to an output directory" bash ./genDiversity_shared.sh
for rg in wholePop Trotter Pacer; do bash OUTPUT_DIR="Path to same output directory" ./genDiversity_per_group.sh "$rg"; done
OUTPUT_DIR="Path to same output directory" bash ./genDiversity_aggregate.sh
```

## Validating outputs

After the pipeline finishes, check that every expected per-group reference file is present with the right header and row count:

```bash
bash scripts/validate_popRefs.sh --mode upstream --root "$OUTPUT_DIR"
```

The same script runs in `--mode popFiles` against a populated GPA `popFiles/` directory (see [`MIGRATION.md`](./MIGRATION.md) §5 for how GPA copies and invokes it).

## Benchmarking the refactor

`scripts/benchmark/run_benchmark.sh` builds a 25-Trotter + 25-Pacer subset from an existing full-run `OUTPUT_DIR` and runs `per_group.sh × 3 + aggregate.sh` end-to-end, finishing with the validator — a fast regression test when changing pipeline logic:

```bash
bash scripts/benchmark/run_benchmark.sh \
     --source results_<timestamp> \
     --target results_benchmark_25x25
```
