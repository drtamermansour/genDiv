# genDiv

This pipeline makes use of 500+ horses' genotyping data to assess genetic diversity in the population and create population-level reference ranges for future breeding programs.

## Environment

[`environment.yml`](./environment.yml) is the single source of truth for the runtime environment. Build from it:

```bash
mamba env create -f environment.yml   # first time
mamba env update -f environment.yml   # after environment.yml changes
conda activate genDiv
```

Two things that file will not do for you:

- **`rclone` is not in it.** The pipeline runs `module load rclone` (`genDiversity_shared.sh:12`), so rclone comes from the cluster module system with its remotes already configured. Installing a conda copy would shadow it. You need a working `remote_UCDavis_GoogleDr` remote before the first run.
- **Do not regenerate it with `conda env export`.** That records one machine's solved transitive closure and erases the distinction between "we need this" and "something else pulled it in" — which is how `r-dplyr`, `r-tidyr`, and `tqdm` once went undeclared while the pipeline kept working locally. Edit `environment.yml` by hand.

To confirm the spec still covers the code after adding or changing a script:

```bash
python3 scripts/check_env.py     # exits non-zero on any undeclared dependency
```

It scans every `.py`, `.R`, and `.sh` in the repo for imports, `library()` calls, and CLI invocations, maps each to its conda package name, and fails on anything `environment.yml` does not declare. It needs only a bare `python3`, so it runs without the environment activated.

## Running

The pipeline expects sibling repositories next to this one (`../InfiniTier/` and `../Horse_parentage_SNPs/equCab3/`) and a configured rclone remote (`remote_UCDavis_GoogleDr`). See `CLAUDE.md` → "Repository layout expectations" for details.

```bash
git clone git@github.com:drtamermansour/genDiv.git
cd genDiv
conda activate genDiv
bash ./genDiversity.sh
```

The pipeline is split across five bash files at the repo root (`genDiversity.sh` is the wrapper; `genDiversity_common.sh` holds shared CONFIG/helpers; `genDiversity_shared.sh`, `genDiversity_per_group.sh`, and `genDiversity_aggregate.sh` do the work). Each subscript is independently runnable, useful for refreshing just one stage:

```bash
## Use the SAME OUTPUT_DIR for all three stages — later stages read earlier stages' files.
export OUTPUT_DIR=results_<timestamp>

bash ./genDiversity_shared.sh
for rg in wholePop Trotter Pacer; do bash ./genDiversity_per_group.sh "$rg"; done
bash ./genDiversity_aggregate.sh
```

## What a run produces

Everything lands in one self-contained `results_<timestamp>/` directory — downloaded inputs and all outputs together, so a run is reproducible and disposable as a unit. Plan for **a few GB of disk and an overnight wall clock**: the reference run `results_20260421_200003` is 2.4 GB and was started in the evening and collected the next morning. It runs 10-way parallel by default (`nthreads` in `genDiversity_common.sh`).

| Subdirectory | Size (reference run) | Contents |
|---|---:|---|
| `divStats/` | 1.3 GB | The deliverables: ROH calls and consensus regions, heterozygosity / F_ROH / F_SNP references, PCA, KING and IBS relatedness, FST, ROH_common, ROH islands, and `ne/` effective-population-size trajectories |
| `rep_ROHRM/` | 48 MB | ROH-based relationship matrices and GRM-vs-ROHRM comparisons, per group and per ROH cutoff |
| `LD_pruned/`, `filtered/`, `preprocess/`, `dedup/`, `inspect/` | ~310 MB | QC and filtering intermediates, kept so any number can be traced back to the genotypes it came from |
| `SNPdata_iScan_Standardbred/` | 713 MB | The rclone-downloaded raw input |
| `run.log` | — | Full stdout+stderr of every invocation into this directory |

Outputs are also rclone-uploaded to Google Drive as they are produced, so results survive the compute node.

The nine per-group reference files that the downstream GPA report pipeline consumes are specified in [`MIGRATION.md`](./MIGRATION.md) — that document, not this one, is the contract.

## Exploratory analyses

`explore.sh` runs one-off analyses against a finished run. It is not part of the pipeline and never runs automatically:

```bash
OUTPUT_DIR=results_<timestamp> bash explore.sh              # all tasks
OUTPUT_DIR=results_<timestamp> bash explore.sh <task_name>  # one task
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

## Upload outputs
If you would like to re-upload the outputs of any run (either actual or benchmarking) to the folder `outputs` in GDrive, you can run this command:

```bash
OUTPUT_DIR=results_<timestamp> bash upload_outputs.sh
``` 
