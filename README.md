# genDiv

This pipeline makes use of 500+ horses' genotyping data to assess genetic diversity in the population and create population-level reference ranges for future breeding programs.

For architecture details (script structure, data flow, parameter definitions) see [`CLAUDE.md`](./CLAUDE.md). Open suggestions and TODOs are tracked in [`REFACTORING.md`](./REFACTORING.md).

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
  bioconda::plink bioconda::plink2 bioconda::bcftools \
  bioconda::beagle bioconda::rclone

## Optional / not currently invoked by genDiversity.sh — keep if you plan to wire them in:
##   bioconda::gcta bioconda::bedtools bioconda::snakemake=9.16.3

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
