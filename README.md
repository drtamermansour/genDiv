# genDiv

This pipeline make use of 500+ horses genotyping data to assess the genetic diversity in the population and create population level reference ranges to be used in future breeding programs


## How to run


```bash
## update the conda package manager itself
conda update conda

## update all packages in an environment
conda update --all

## clean up unused packages and caches
conda clean --all

## create a new conda environment for genomic related GWAS tools
mamba create -n genDiv conda-forge::r-base=4.5.2 conda-forge::r-ggplot2=4.0.1 conda-forge::r-gridextra=2.3 \
              conda-forge::r-qqman=0.1.9 conda-forge::r-viridis=0.6.5 conda-forge::r-reshape2=1.4.5 \
              conda-forge::r-ggally=2.4.0 conda-forge::openpyxl conda-forge::scikit-allel conda-forge::r-effsize=0.8.1 \
              conda-forge::numpy conda-forge::pandas \
              bioconda::plink bioconda::plink2 bioconda::bcftools bioconda::gcta bioconda::bedtools bioconda::beagle bioconda::snakemake=9.16.3

## Added 
## conda-forge::numpy conda-forge::scikit-allel \
## Missing        
## conda-forge r-hierfstat=0.5_11 
## conda-forge matplotlib=3.10.8 seaborn=0.13.2 scipy=1.17.0
conda activate genDiv


## Create the working directory of the project
git clone git@github.com:drtamermansour/genDiv.git && cd genDiv
bash ./genDiversity.sh
```