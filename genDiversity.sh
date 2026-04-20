#!/usr/bin/env bash
set -eo pipefail

# Load CONFIG, helper functions, working-dir vars, and log redirection.
source "$(dirname "$0")/genDiversity_common.sh"


## Run whole-pop preprocessing (Sections 1–3).
bash "$(dirname "$0")/genDiversity_shared.sh"

## Re-declare canonical paths produced by shared.sh (it runs in a subshell, so
## its variables don't carry over). All of these are deterministic functions of
## $OUTPUT_DIR.
SNPdata="$(pwd)/${OUTPUT_DIR}/SNPdata_iScan_Standardbred"
docs="$(pwd)/${OUTPUT_DIR}/Miscellaneous_documents_standardbred"
pl1_filtered="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered"
vcf_filtered="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.vcf"
pl1_pruned="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune"
vcf_pruned="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.LD_prune.vcf"

############## Stats on diversity ##################
log "Section 5: Diversity statistics"
mkdir -p ${OUTPUT_DIR}/divStats

##########################################
## PCA Assessment
##########################################
## PCAs are called "loadings" because they represent the weights or coefficients that determine how much each original variable "loads" onto or contributes to a specific PC.
pca_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.pca"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --real-ref-alleles --autosome --pca 'allele-wts' \
       --output-chr 'chrM' --out "$pca_prefix"

rclone -v copy ${OUTPUT_DIR}/divStats --include "filtered.LD_prune.pca.eigen*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

Rscript -e 'args=(commandArgs(TRUE));'\
'val <- read.table(paste(args[1],"eigenval",sep="."));'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = args[2]);'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "principal Component", ylab = "Variance explained in %");'\
'dev.off();' "$pca_prefix" "${OUTPUT_DIR}/divStats/Var_PCs.jpg"
rclone -v copy ${OUTPUT_DIR}/divStats/Var_PCs.jpg "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Color samples on the PCA plots by Sex
awk 'BEGIN{FS=OFS="\t";a["IID"]="sex"}NR==FNR{if($5==1)a[$2]="male";else a[$2]="female";next}{print $0,a[$2]}' $pl1_pruned.fam $pca_prefix.eigenvec > $pca_prefix.eigenvec.wSex
eigenvec_suffix="wSex"; color_column="sex"; out_png="${OUTPUT_DIR}/divStats/pca_plot_sex.png";
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" 6 factor
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Color samples on the PCA plots by Gait
awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait $pca_prefix.eigenvec > $pca_prefix.eigenvec.wGait
eigenvec_suffix="wGait"; color_column="Gait"; out_png="${OUTPUT_DIR}/divStats/pca_plot_Gait.png";
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" 6 factor
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Color samples on the PCA plots by Book Size
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize $pca_prefix.eigenvec > $pca_prefix.eigenvec.wBook_Size
eigenvec_suffix="wBook_Size"; color_column="Book_Size"; out_png="${OUTPUT_DIR}/divStats/pca_plot_BookSize.png";
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" 6 factor
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me


## Identify Trotter samples segregating on PC2
## After exlcusion of highly related animal, this subpopulation is segregating on PC4 (I keep the name of file on_PC2 to avoid confusion ) 
cat $pca_prefix.eigenvec | awk 'BEGIN{FS=OFS="\t"}{if($6>0.1)print $2}' | grep -Fwf - $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv > ${OUTPUT_DIR}/divStats/Trotters_segregating_on_PC2.csv || true
rclone -v copy ${OUTPUT_DIR}/divStats/Trotters_segregating_on_PC2.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
## Identify Pacer samples co-segregating with Trotters on PC1
## 10 samples; currently labeled as undefined.  
cat $pca_prefix.eigenvec | awk 'BEGIN{FS=OFS="\t"}{if($3<0)print $2}' | grep -Fwf - $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv | grep "Pacer" > ${OUTPUT_DIR}/divStats/Pacers_cosegregating_withTrotters_on_PC1.csv || true
rclone -v copy ${OUTPUT_DIR}/divStats/Pacers_cosegregating_withTrotters_on_PC1.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
## Identify Trotter samples co-segregating with Pacers on PC1
## 8 samples; currently labeled as undefined.  
cat $pca_prefix.eigenvec | awk 'BEGIN{FS=OFS="\t"}{if($3>0)print $2}' | grep -Fwf - $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv | grep "Trotter" > ${OUTPUT_DIR}/divStats/Trotters_cosegregating_withPacers_on_PC1.csv || true
rclone -v copy ${OUTPUT_DIR}/divStats/Trotters_cosegregating_withPacers_on_PC1.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me


## 2A. PCA (Trotter only)
pca_prefix_trot="${OUTPUT_DIR}/divStats/filtered.LD_prune.Trotter.pca"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --keep <(grep "Trotter" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait) \
       --autosome --pca \
       --output-chr 'chrM' --out "$pca_prefix_trot"

rclone -v copy ${OUTPUT_DIR}/divStats --include "filtered.LD_prune.Trotter.pca.eigen*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

Rscript -e 'args=(commandArgs(TRUE));'\
'val <- read.table(paste(args[1],"eigenval",sep="."));'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = args[2]);'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "principal Component", ylab = "Variance explained in %");'\
'dev.off();' "$pca_prefix_trot" "${OUTPUT_DIR}/divStats/Var_PCs.Trotter.jpg"
rclone -v copy ${OUTPUT_DIR}/divStats/Var_PCs.Trotter.jpg "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Color samples on the PCA plots by Book Size
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize $pca_prefix_trot.eigenvec > $pca_prefix_trot.eigenvec.wBook_Size
eigenvec_suffix="wBook_Size"; color_column="Book_Size"; out_png="${OUTPUT_DIR}/divStats/pca_plot_BookSize.Trotter.png";
Rscript scripts/pca_plots.R "$pca_prefix_trot" "$eigenvec_suffix" "$color_column" "$out_png" 3 factor
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me


## 2B. PCA (Pacer only)
pca_prefix_pace="${OUTPUT_DIR}/divStats/filtered.LD_prune.Pacer.pca"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --keep <(grep "Pacer" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait) \
       --autosome --pca \
       --output-chr 'chrM' --out "$pca_prefix_pace"

rclone -v copy ${OUTPUT_DIR}/divStats --include "filtered.LD_prune.Pacer.pca.eigen*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

Rscript -e 'args=(commandArgs(TRUE));'\
'val <- read.table(paste(args[1],"eigenval",sep="."));'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = args[2]);'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "principal Component", ylab = "Variance explained in %");'\
'dev.off();' "$pca_prefix_pace" "${OUTPUT_DIR}/divStats/Var_PCs.Pacer.jpg"
rclone -v copy ${OUTPUT_DIR}/divStats/Var_PCs.Pacer.jpg "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Color samples on the PCA plots by Book Size
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize $pca_prefix_pace.eigenvec > $pca_prefix_pace.eigenvec.wBook_Size
eigenvec_suffix="wBook_Size"; color_column="Book_Size"; out_png="${OUTPUT_DIR}/divStats/pca_plot_BookSize.Pacer.png";
Rscript scripts/pca_plots.R "$pca_prefix_pace" "$eigenvec_suffix" "$color_column" "$out_png" 3 factor
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

########################################################
## 1. Effective number of alleles (\(A_{e}\)) 
## A_e represents the number of equally frequent alleles required to achieve the same level of expected heterozygosity (\(H_{e}\)) observed in a population
########################################################

## Formula
## \(A_{e} = \frac{1}{\sum p_{i}^{2}}\)
## where \(p_{i}\) is the frequency of the \(i^{th}\) allele

## Example Calculation
## Suppose a single locus has three alleles with the observed frequencies (0.6, 0.3, 0.1) in a population.
## 1. Calculate the squared frequencies: 0.36, 0.09, and 0.01
## 2. Calculate \(A_{e}\): 1/(0.36 + 0.09 + 0.01) = 1/0.46 = 2.17
## This result means that although there are 3 distinct alleles, the population's genetic diversity is equivalent to a population with only 2.17 equally frequent alleles. 

## Calculate \(A_{e}\) for each SNP
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --freq \
    --out "$pl1_pruned.freq_stats"
awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.afreq" > "$pl1_pruned.freq_stats.wholePop.afreq.Ae"
awk -v pop="wholePop" 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
    { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"pop, mean_Ae, "SD_Ae_in_"pop, sd_Ae } }' "$pl1_pruned.freq_stats.wholePop.afreq.Ae"
#Mean_Ae_in_wholePop     1.5616  SD_Ae_in_wholePop       0.321994

## Calculate \(A_{e}\) for each SNP per gait subpopulation
group="gait"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group \
    --loop-cats 'PHENO1' --freq \
    --out "$pl1_pruned.freq_stats"
#--loop-cats: Processing category 'Pacer' (271 samples).
#--loop-cats: Processing category 'Trotter' (271 samples).

## Calculate Mean \(A_{e}\) and standard deviation per gait subpopulation
for pop in Trotter Pacer; do
    awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.$pop.afreq" > "$pl1_pruned.freq_stats.$pop.afreq.Ae"
    awk -v pop=$pop 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
        { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"pop, mean_Ae, "SD_Ae_in_"pop, sd_Ae } }' "$pl1_pruned.freq_stats.$pop.afreq.Ae"
done
#Mean_Ae_in_Trotter      1.53063 SD_Ae_in_Trotter        0.339752
#Mean_Ae_in_Pacer        1.55017 SD_Ae_in_Pacer  0.331251

## Calculate \(A_{e}\) for each SNP per book size in each gait subpopulation
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize \
    --loop-cats 'PHENO1' --freq \
    --out "$pl1_pruned.freq_stats"
#--loop-cats: Processing category 'Pacer_HIGH' (75 samples).
#--loop-cats: Processing category 'Pacer_LOW' (90 samples).
#--loop-cats: Processing category 'Pacer_MEDIUM' (106 samples).
#--loop-cats: Processing category 'Trotter_HIGH' (58 samples).
#--loop-cats: Processing category 'Trotter_LOW' (96 samples).
#--loop-cats: Processing category 'Trotter_MEDIUM' (117 samples).

## Calculate Mean \(A_{e}\) and standard deviation per book size
for pop in Trotter Pacer; do
    for book in LOW MEDIUM HIGH;do
        awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.${pop}_${book}.afreq" > "$pl1_pruned.freq_stats.${pop}_${book}.afreq.Ae"
        awk -v gp=${pop}_${book} 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
            { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"gp, mean_Ae, "SD_Ae_in_"gp, sd_Ae } }' "$pl1_pruned.freq_stats.${pop}_${book}.afreq.Ae"
    done
done
#Mean_Ae_in_Pacer_LOW    1.55141 SD_Ae_in_Pacer_LOW      0.330248
#Mean_Ae_in_Pacer_MEDIUM 1.54804 SD_Ae_in_Pacer_MEDIUM   0.331937
#Mean_Ae_in_Pacer_HIGH   1.53805 SD_Ae_in_Pacer_HIGH     0.338842
#Mean_Ae_in_Trotter_LOW  1.53338 SD_Ae_in_Trotter_LOW    0.338901
#Mean_Ae_in_Trotter_MEDIUM       1.53176 SD_Ae_in_Trotter_MEDIUM 0.340108
#Mean_Ae_in_Trotter_HIGH 1.50661 SD_Ae_in_Trotter_HIGH   0.34827


Rscript scripts/effAllele_stats.R "${OUTPUT_DIR}" &> ${OUTPUT_DIR}/divStats/effAllele_stats.txt
Rscript scripts/plot_Ae.R "${OUTPUT_DIR}"
rclone -v copy ${OUTPUT_DIR}/divStats/effAllele_stats.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Ae/" --drive-shared-with-me
rclone -v copy ${OUTPUT_DIR}/divStats/Figure_Ae_BookSize.tiff "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Ae/" --drive-shared-with-me
#rclone -v copy ${OUTPUT_DIR}/divStats/Figure_Ae_BookSize.pdf "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Ae/" --drive-shared-with-me

##########################################
## 2. Fst between subpopulations (genders, gait types, and book sizes)
##########################################
## The fixation index can range from 0 to 1, where 0 means complete sharing of genetic material and 1 means no sharing. 
## For values equal to 1(meaning no sharing), scientists say that the populations are fixed.
## Effects of marker type and filtering criteria on QST-FST comparisons: https://pmc.ncbi.nlm.nih.gov/articles/PMC6894560/
for group in sex gait bookSize; do
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group \
        --fst 'PHENO1' 'blocksize=2000' \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_$group
done

find ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_*.summary -maxdepth 1 -type f | grep -v "\.x\." | xargs cat > ${OUTPUT_DIR}/divStats/autosomal.fst.summary
rclone -v copy ${OUTPUT_DIR}/divStats/autosomal.fst.summary "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me

group="bookSize"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep <(grep "Trotter" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait) \
    --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group \
    --fst 'PHENO1' 'blocksize=2000' \
    --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_$group.Trotter
rclone -v copy ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_bookSize.Trotter.fst.summary "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me

plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep <(grep "Pacer" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait) \
    --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group \
    --fst 'PHENO1' 'blocksize=2000' \
    --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_$group.Pacer
rclone -v copy ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_bookSize.Pacer.fst.summary "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me

Rscript scripts/fst_stats.R "${OUTPUT_DIR}" &> ${OUTPUT_DIR}/divStats/fst_stats.txt
rclone -v copy ${OUTPUT_DIR}/divStats/fst_stats.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me

##########################################
## 3. Expected and observed heterozygosity and inbreeding coefficient
##########################################
## An inbreeding coefficient (COI) is a measure of the probability that an individual will have two copies of an allele that are identical by descent from a common ancestor.
## A higher COI means more predictability of traits but also a greater risk of genetic health problems due to inbreeding depression
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --het 'cols=fid,hom,het,nobs,f' \
    --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                END { for(i=bmin;i<=bmax;++i){if(i==0) print -1*size,size,a[i]/1;else if(i<0) print (i-1)*size,i*size,a[i]/1;else print i*size,(i+1)*size,a[i]/1 }}'  <(tail -n+2 ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het) > ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.histo

## generate a summary table of heterozygosity and inbreeding coefficient in the two subpopulations and the whole cohort
awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' \
     ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het > ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait

INPUT_HET="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait.sumStats.csv"
python scripts/summary_het.py -i "$INPUT_HET" -o "$OUTPUT_FILE"

## generate a summary table of heterozygosity and inbreeding coefficient in the three book size in the two subpopulations and the whole cohort
awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' \
     ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het > ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize

INPUT_HET="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize.sumStats.csv"
python scripts/summary_het.py -i "$INPUT_HET" -o "$OUTPUT_FILE"

rclone -v copy ${OUTPUT_DIR}/divStats --include "filtered.LD_prune.het_stats.het*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/het_and_COI/" --drive-shared-with-me

## PCA: Color samples on the PCA plots by COI
# Whole population
pca_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.pca"
awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het) $pca_prefix.eigenvec > $pca_prefix.eigenvec.wCOI
eigenvec_suffix="wCOI"; color_column="COI"; out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.png";
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" 6 numeric
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Trotters only
pca_prefix_trot="${OUTPUT_DIR}/divStats/filtered.LD_prune.Trotter.pca"
awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het) $pca_prefix_trot.eigenvec > $pca_prefix_trot.eigenvec.wCOI
eigenvec_suffix="wCOI"; color_column="COI"; out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.Trotter.png";
Rscript scripts/pca_plots.R "$pca_prefix_trot" "$eigenvec_suffix" "$color_column" "$out_png" 3 numeric
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Pacers only
pca_prefix_pace="${OUTPUT_DIR}/divStats/filtered.LD_prune.Pacer.pca"
awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het) $pca_prefix_pace.eigenvec > $pca_prefix_pace.eigenvec.wCOI
eigenvec_suffix="wCOI"; color_column="COI"; out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.Pacer.png";
Rscript scripts/pca_plots.R "$pca_prefix_pace" "$eigenvec_suffix" "$color_column" "$out_png" 3 numeric
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me


##########################################
# 4. Runs of homozygosity (ROH)
##########################################
group="gait"
case="Trotter"

##########################################
## 4A. ROH using Plink (Pruned dataset) -- This is not actual useful. Only for testing the effect of pruning and using fewer count of markers 
##########################################
## Plink --homozyg
## By default, only runs of homozygosity containing at least 100 SNPs, and of total length ≥ 1000 kilobases, are noted. You can change these minimums with --homozyg-snp and --homozyg-kb, respectively.
## By default, a ROH must have at least one SNP per 50 kb on average; change this bound with --homozyg-density.
## By default, if two consecutive SNPs are more than 1000 kb apart, they cannot be in the same ROH; change this bound with --homozyg-gap.
## By default, a ROH can contain an unlimited number of heterozygous calls; you can impose a limit with --homozyg-het. (This flag was silently ignored by PLINK 1.07.)
## By default, the scanning window contains 50 SNPs; change this with --homozyg-window-snp.
## By default, a scanning window hit can contain at most 1 heterozygous call and 5 missing calls; change these limits with --homozyg-window-het and --homozyg-window-missing, respectively.
## By default, for a SNP to be eligible for inclusion in a ROH, the hit rate of all scanning windows containing the SNP must be at least 0.05; change this threshold with --homozyg-window-threshold.
## Due to how the scanning algorithm works, it is possible for a reported run of homozygosity to be adjacent to a few unincluded homozygous variants. This is generally harmless, but if you wish to extend the ROH to include them, use the 'extend' modifier. (Note that the --homozyg-density bound can prevent extension, and --homozyg-gap affects which variants are considered adjacent.)

plink --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group $case \
        --homozyg 'extend' \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group


## Outputs of --homozyg include: 
## .hom (run-of-homozygosity list)
##  FID, IID, PHE (henotype value), CHR, SNP1, SNP2, POS1, POS2, KB, NSNP, DENSITY, PHOM (% of homozygous), PHET (% of heterozygous), PMISS (% of missing), SENSITIVITY, SPECIFICITY

## .hom.indiv (sample-based runs-of-homozygosity report) which has the following columns: 
## FID	Family ID
## IID	Within-family ID
## PHE	Phenotype value
## NSEG	Number of runs of homozygosity
## KB	Total length of runs (kb)
## KBAVG	Average length of runs (kb)
## Calculate summary stats from .hom.indiv 
awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/count, sum5/count, sum6/count }' ${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 15.63 
##Average of the total length of runs (kb) across all samples: 164,427.23
##Average of the average length of runs (KBAVG) across all samples: 10,439.76


## Rscript that plots the correlation between "KB" and "KBAVG" from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group.hom"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_pruned/" --drive-shared-with-me

##########################################
## 4B. ROH using Plink (Filtered dataset without LD pruning)
##########################################
plink --bfile "$pl1_filtered" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group $case \
        --homozyg 'extend' \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group

awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/count, sum5/count, sum6/count }' ${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 33.35
##Average of the total length of runs (kb) across all samples: 345,017.70
##Average of the average length of runs (KBAVG) across all samples: 10,305.44


## Rscript that plots the correlation between  KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group.hom"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_${OUTPUT_DIR}/filtered/" --drive-shared-with-me

##########################################
## 4C. ROH using Plink (Filtered dataset without LD pruning (with group option))
##########################################
plink --bfile "$pl1_filtered" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group $case \
        --homozyg group 'extend' \
        --homozyg-window-snp 20 \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group

awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/count, sum5/count, sum6/count }' ${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 35.46
##Average of the total length of runs (kb) across all samples: 364,115.99
##Average of the average length of runs (KBAVG) across all samples: 10,211.80

## Rscript that plots the correlation between  KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group.hom"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_filtered_gp/" --drive-shared-with-me

############################################
## 4D. Studying ROH using Howard et al. (2016) approach -- This was for testing only and is not currently in use
############################################
mkdir -p $work_dir/${OUTPUT_DIR}/rep_ROHRM
rohrm_dir="$work_dir/${OUTPUT_DIR}/rep_ROHRM"

## ROH Analysis with Sub-populations and Phenotypes
## Having a headerless tab-separated file with 3 columns: The 2nd column has subject ids matching the VCF and the 3rd column has the a binary phenotype, let us do the following:
## 1. ROH Calling (Per Individual): We will implement a scanner that checks hap1 == hap2. Any contiguous stretch of matching haplotypes longer than the cutoff (e.g., 1 Mb) is flagged as an ROH.
## 2. Island Detection: We will calculate the frequency of ROHs at every SNP, find the "Top 5%" cutoff, and merge contiguous high-frequency SNPs into "Islands".
##    i.e., an "ROH Island" is defined strictly as a contiguous block of SNPs where every single SNP is in the Top 5% of frequencies.
## 3. Phenotype Integration: the analysis will be repeated for sub-populations defined in the phenotype file).
## 4. make a plot to show the ROH frequency across the genome for the whole population and each sub-population.
## 5. Calc the average (±SD) proportion of the genome in a ROH for the whole population and each sub-population.

roh_mb_cutoff=1.0  # in Megabases (Mb)
phenotypes="${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait"
python $scripts/ROH_analysis.py $vcf_filtered.norm.phased.vcf.gz $phenotypes $roh_mb_cutoff "$rohrm_dir" > $rohrm_dir/ROH_analysis.$roh_mb_cutoff.log

rclone -v copy $rohrm_dir/ROH_Frequency_Plot.$roh_mb_cutoff.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
rclone -v copy $rohrm_dir/ROH_Islands_Detailed.$roh_mb_cutoff.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
#rclone -v copy $rohrm_dir/ROH_Subpop_Stats.$roh_mb_cutoff.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
rclone -v copy $rohrm_dir/ROH_analysis.$roh_mb_cutoff.log "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me


tail -n+2 $rohrm_dir/ROH_Islands_Detailed.1.0.csv | awk -F"," '{sum += $5; array[NR] = $5} END {
    mean = sum / NR;
    for (x=1; x<=NR; x++) {
        sumsq += ((array[x] - mean)**2);
    }
    std_dev = sqrt(sumsq / NR); # Population standard deviation
    # For sample standard deviation, use sqrt(sumsq / (NR - 1)) if NR > 1
    min_snp_cutoff = mean - (2 * std_dev)

    print "Mean: " mean;
    print "Standard Deviation: " std_dev
    print "2 SD below Mean: " min_snp_cutoff
}' ## Mean: 28.6454 // Standard Deviation: 41.1817 // 2 SD below Mean: -53.7179 //There is no need to filter based on this criterion as it results in a negative value. As an alternative, we can use 3 as a minimum SNP count threshold for defining ROH islands.

## Filter ROH islands with at least 3 SNPs
awk -F"," 'NR==1 || $5 >= 3' $rohrm_dir/ROH_Islands_Detailed.1.0.csv > $rohrm_dir/Filtered_ROH_Islands_Detailed.1.0.csv
## for each sub-population in column 1, calculate the number of ROH regions (NR), the total length of ROH islands (sum of $4-$3), the average length, and the maximum and average SNP frequency (column 6)
awk -F"," 'NR>1 && $1!="" {
    grp=$1
    if (!(grp in seen)) { seen[grp]=1; order[++norder]=grp }
    len = ($4 - $3) + 0
    freq = ($6 + 0)
    count[grp]++
    sumlen[grp] += len
    sumfreq[grp] += freq
    if (!(grp in maxfreq) || freq > maxfreq[grp]) maxfreq[grp] = freq
}
END {
    print "Group,NR,Total_len_bp,Avg_len_bp,Max_freq,Avg_freq"
    for(i=1;i<=norder;i++) {
        g = order[i]
        printf "%s,%d,%.0f,%.2f,%.6f,%.6f\n", g, count[g], sumlen[g], sumlen[g]/count[g], maxfreq[g], sumfreq[g]/count[g]
    }
}' $rohrm_dir/Filtered_ROH_Islands_Detailed.1.0.csv > $rohrm_dir/Filtered_ROH_Subpop_Stats.csv

rclone -v copy $rohrm_dir/Filtered_ROH_Subpop_Stats.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
#python $scripts/ROH_analysis_withFilter.py filtered.norm.phased.vcf.gz phenotypes.txt $roh_mb_cutoff > Filtered_ROH_analysis.$roh_mb_cutoff.log

##########################################
# 4E. ROH using bcftools/roh (Filtered dataset without LD pruning) -- This is the final approved approach
##########################################
## Run bcftools roh
bcftools roh -G30 --estimate-AF - $vcf_filtered.norm.phased.vcf.gz -o ${OUTPUT_DIR}/divStats/roh_out.txt
##Number of target samples: 560
##Number of --estimate-AF samples: 560
##Number of sites in the buffer/overlap: unlimited
##Number of lines total/processed: 57829/57829 (old: 58106/58106)
##Number of lines ${OUTPUT_DIR}/filtered/no AF/no alt/multiallelic/dup: 0/0/0/0/0

grep -E "^RG|^#" ${OUTPUT_DIR}/divStats/roh_out.txt > ${OUTPUT_DIR}/divStats/roh_out_RG.txt
## Summary stats by RG
awk 'BEGIN{print "IID\tNSEG\tKB\tKBAVG"} $1=="RG"{n[$2]++; sum[$2]+=$6} END{for (s in n) printf "%s\t%d\t%.2f\t%.2f\n", s, n[s], sum[s]/1000, (sum[s]/1000)/n[s]}' ${OUTPUT_DIR}/divStats/roh_out_RG.txt > ${OUTPUT_DIR}/divStats/roh_summary_by_RG.txt
awk 'NR > 1{ sum2 += $2; sum3 += $3; sum4 += $4 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum2/count, sum3/count, sum4/count }' ${OUTPUT_DIR}/divStats/roh_summary_by_RG.txt
##Average Number of runs of homozygosity (NSEG) : 86.53
##Average of the total length of runs (kb) across all samples: 456,016.51
##Average of the average length of runs (KBAVG) across all samples: 5,256.67

## filtration to match the PLINK quality suggestions 
#Minimum ROH length (--homozyg-kb) 1000 kb
awk '/^#/ || $6 >= 1000000' ${OUTPUT_DIR}/divStats/roh_out_RG.txt > ${OUTPUT_DIR}/divStats/roh.L1.txt
#Minimum number of SNPs in ROH (--homozyg-snp) 50
awk '/^#/ || $7 >= 50' ${OUTPUT_DIR}/divStats/roh.L1.txt > ${OUTPUT_DIR}/divStats/roh.L2.txt
#Quality scores
awk -v size=2 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(grep -v "^#" ${OUTPUT_DIR}/divStats/roh.L2.txt) > ${OUTPUT_DIR}/divStats/roh.L2.histo 
awk '/^#/ || $8 >= 20' ${OUTPUT_DIR}/divStats/roh.L2.txt > ${OUTPUT_DIR}/divStats/roh.L3.txt

## Summary stats by RG after QC filtration && Stratify the file by the gait type
## The stats will be recalculated again later with Froh using "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt"
awk 'BEGIN{print "IID\tNSEG\tKB\tKBAVG"} $1=="RG"{n[$2]++; sum[$2]+=$6} END{for (s in n) printf "%s\t%d\t%.2f\t%.2f\n", s, n[s], sum[s]/1000, (sum[s]/1000)/n[s]}' ${OUTPUT_DIR}/divStats/roh.L3.txt > ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt
awk 'BEGIN{FS=OFS="\t";gait["IID"]="gait"}FNR==NR{gait[$2]=$3;next} {if(gait[$1])print $0,gait[$1];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt > ${OUTPUT_DIR}/divStats/roh.L3_gait.txt
INPUT_ROH="${OUTPUT_DIR}/divStats/roh.L3_gait.txt"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/roh.L3_gait.sumStats.csv"
python scripts/summary_roh.py -i "$INPUT_ROH" -o "$OUTPUT_FILE"
rclone -v copy ${OUTPUT_DIR}/divStats/roh.L3_gait.sumStats.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
#Subgroup,          N,      NSEG,           KB,                         KBAVG
#Whole Population,  560,    54.58 +/- 9.43, 413086.61 +/- 103171.97,    7533.22 +/- 1210.63
#Pacer,             271,    50.34 +/- 7.11, 379860.14 +/- 81886.48,     7535.88 +/- 1224.94
#Trotter,           271,    59.28 +/- 9.02, 451519.43 +/- 105452.45,    7580.75 +/- 1166.36
#undefined,         18,     47.67 +/- 11.92,334702.03 +/- 138748.03,    6777.52 +/- 1454.36


## Rscript that plots the correlation between KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
## Similar analysis will be done later after calculation of related matrices
roh_indiv="${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_summary_by_RG_L3"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

## A per-base consensus ROH where ≥25% of samples are in ROH filtered by minimum size 500 kb and stratified by gait type
## With and without applying a smoothing function to the per-base coverage data to reduce noise before identifying consensus ROH regions
#roh_RG=${OUTPUT_DIR}/divStats/roh_out_RG
roh_RG=${OUTPUT_DIR}/divStats/roh.L3
# 1. Convert RG output → BED format
awk 'BEGIN{OFS="\t"} $1=="RG" {print $3, $4-1, $5, $2}' "${roh_RG}.txt" > "${roh_RG}.bed"
# 2. Ensure ROHs from the same sample do not double-count
cut -f4 "${roh_RG}.bed" | sort -u | while read S; do
  awk -v s="$S" '$4==s' "${roh_RG}.bed" | sort -k1,1 -k2,2n | bedtools merge -i - | awk -v s="$S" 'BEGIN{OFS="\t"}{print $1,$2,$3,s}'
done | sort -k1,1 -k2,2n > "${roh_RG}.merged_per_sample.wholePop.bed"

# subset the bed file for each subpopulation
for rg in "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do 
    grep "$rg" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize | cut -f2 | grep -f - "${roh_RG}.merged_per_sample.wholePop.bed" > "${roh_RG}.merged_per_sample.${rg}.bed"
done

# 3. Calculate per-base ROH frequency (i.e., how many samples are in ROH at each base position)
awk '$1 ~ /^[0-9]+$/' $reference_fai | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2}' > ${OUTPUT_DIR}/divStats/autosomes.genome
for rg in "wholePop" "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do 
    bedtools genomecov -i "${roh_RG}.merged_per_sample.${rg}.bed" -g ${OUTPUT_DIR}/divStats/autosomes.genome -bg > "${roh_RG}.per_base_coverage.${rg}.bed"
    awk -v size=5 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($4/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                      END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  "${roh_RG}.per_base_coverage.${rg}.bed" > "${roh_RG}.per_base_coverage.${rg}.histo"
    # upload bed files
    # pause for now to save space
    #rclone -v copy ${roh_RG}.per_base_coverage.${rg}.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/freq/" --drive-shared-with-me

    # upload histo files
    rclone -v copy ${roh_RG}.per_base_coverage.${rg}.histo "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/freq/" --drive-shared-with-me
done

# 4. Identify consensus ROH regions (≥25% of samples in ROH) and merge adjacent regions (minimum size 500 kb)
# With and without appling a smoothing function which adjust the per-base coverage value of regions briding intervals with high coverage. The function would assign the average coverage of the region and two flanking regions to the bridged interval.
for rg in "wholePop" "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do
  num_samples=$(cut -f4 "${roh_RG}.merged_per_sample.${rg}.bed" | sort -u | wc -l)
  threshold=$(echo "$pct * $num_samples / 100" | bc -l)
  ## Find consensus before smoothing
  awk -v threshold=$threshold 'BEGIN{OFS="\t"} $4 >= threshold {print}' "${roh_RG}.per_base_coverage.${rg}.bed" > "${roh_RG}.consensus_${pct}pct.${rg}.bed"
  bedtools merge -i "${roh_RG}.consensus_${pct}pct.${rg}.bed" -c 4 -o mean | awk -v min_mb=$CONSENSUS_MIN_MB 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=min_mb)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${rg}.bed"

  # Summary stats of consensus ROH regions
  echo "==== consensus ROH in ≥${pct}% of ${rg} samples BEFORE SMOOTHING ======"
  awk -v rg="$rg" -v nsam="$num_samples" 'BEGIN{OFS=",";maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
    {print rg,"\nNo. of segments","Total length (KB)","Ave. length (KB)","Max % of samples in consensus","Average % of samples in consensus",\
    "\n"NR,sumLen,sumLen/NR,(maxConsen/nsam)*100"%",((sumSamples/NR)/nsam)*100"%"}' "${roh_RG}.consensus_${pct}pct.merged.${rg}.bed"

  # Smooth per-base coverage: only average a middle interval if it exactly bridges two adjacent intervals
  # and both flanking intervals are >= threshold while the middle < threshold.
  awk -v thr="$threshold" 'BEGIN{OFS="\t"} {chr[NR]=$1; st[NR]=$2; en[NR]=$3; cov[NR]=$4} END{
      for(i=1;i<=NR;i++){
          newcov=cov[i]
          if(i>1 && i<NR){
              # check perfect contiguity: prev_end == cur_start && cur_end == next_start
              if(en[i-1]==st[i] && en[i]==st[i+1]){
                  if(cov[i-1] >= thr && cov[i+1] >= thr && cov[i] < thr){
                      newcov = (cov[i-1] + cov[i] + cov[i+1]) / 3
                  }
              }
          }
          printf "%s\t%d\t%d\t%.6f\n", chr[i], st[i], en[i], newcov
      }
  }' "${roh_RG}.per_base_coverage.${rg}.bed" | awk 'BEGIN{OFS="\t"}{$4=$4+0;print}' > "${roh_RG}.per_base_coverage.${rg}.smoothed.bed"
  # Generate new consensus using the smoothed per-base coverage
  awk -v threshold="$threshold" 'BEGIN{OFS="\t"} $4 >= threshold {print}' "${roh_RG}.per_base_coverage.${rg}.smoothed.bed" > "${roh_RG}.consensus_${pct}pct.${rg}.smoothed.bed" ## recovered 72 more regions for wholePop
  bedtools merge -i "${roh_RG}.consensus_${pct}pct.${rg}.smoothed.bed" -c 4 -o mean | awk -v min_mb=$CONSENSUS_MIN_MB 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=min_mb)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed"
  
  # upload bed files
  # pause for now to save space
  #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${rg}.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
  #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

  # Summary stats of consensus ROH regions
  echo "==== consensus smoothed ROH in ≥${pct}% of ${rg} samples AFTER SMOOTHING ======"
  awk -v rg="$rg" -v nsam="$num_samples" 'BEGIN{OFS=",";maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
   {print rg,"\nNo. of segments","Total length (KB)","Ave. length (KB)","Max % of samples in consensus","Average % of samples in consensus",\
    "\n"NR,sumLen,sumLen/NR,(maxConsen/nsam)*100"%",((sumSamples/NR)/nsam)*100"%"}' "${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed"
 
done > ${roh_RG}.consensus_${pct}pct.summary.txt
rclone -v copy ${roh_RG}.consensus_${pct}pct.summary.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
grep -A3 "AFTER SMOOTHING" ${roh_RG}.consensus_${pct}pct.summary.txt

## intersect ROH regions of each sample aganist the consensus ROH regions (in each "rg")
## Output: the census length and percentage in each sample (file for each "rg")
roh_RG="${OUTPUT_DIR}/divStats/roh.L3"
for rg in "wholePop" "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do
    consensus_bed=${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed
    consensus_size=$(awk 'BEGIN{sum=0} {sum+=($3-$2)} END {print sum}' ${consensus_bed})
    bed_perSample="${roh_RG}.merged_per_sample.${rg}.bed"   ## no need to loop on ${rg} here. It should be the same if you always used "wholePop" 
    echo -e "IID\tTotal_ROH_in_Consensus_region(bp)\tPercent_of_Consensus_ROH" > ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt
    cut -f4 "${bed_perSample}" | sort -u | while read S; do
      awk -v s="$S" '$4==s' "${bed_perSample}" | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b "${consensus_bed}" | awk -v s="$S" -v cs="$consensus_size" 'BEGIN{OFS="\t"}{size+=($3-$2)} END {print s, size, (size/cs)*100}'
    done >> ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt
    #rclone -v copy ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
done

## merge the intersection (ROH_share) with Trotter/Pacer consensus
## ouput: the census length and percentage in each sample aganist its own gait consensus
head -n1  ${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt > ${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt
for rg in "Trotter" "Pacer";do 
    tail -n+2 ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt;
done >> ${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt  

## merge the intersection (ROH_share) with Trotter_booksize/Pacer_booksize consensus
## ouput: the census length and percentage in each sample aganist its own gait_bookSize consensus
head -n1  ${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt > ${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt
for rg in "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH";do 
    tail -n+2 ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt;
done >> ${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt  

############################################
## 5. F_ROH statistic (currently calculated based on bacftools roh)
############################################
## F_ROH is an inbreeding coefficient based on runs of homozygosity
## Standard practice is to calculate F_{ROH} statistics on all valid ROHs (>1Mb), while restricting "Islands" (signatures of selection) to only the most robust regions.
## per-sample F_ROH = (sum length of ROH for that individual) / (total autosomal genome length).

## calculate the effective autosomal genome length
awk '{print $1"\t"$4}' "$pl1_filtered".bim | grep "^chr" | grep -v "^chrX" > "$pl1_filtered".snp_pos.txt
aut_len=$(sort -k1,1 -k2,2n "$pl1_filtered".snp_pos.txt | \
        awk '{if ($1 == prev_chr) { gap = $2 - prev_pos; \
              if(gap > 0) {if (gap > 1000000) gap = 1000000; total += gap; }}\
              prev_chr=$1; prev_pos=$2} END {print total}') ## 2,261,547,402
echo $aut_len > ${OUTPUT_DIR}/divStats/effective_autosomal_genome_length.txt

awk -v aut_len=$aut_len 'BEGIN{FS=OFS="\t";}NR==1{print $0,"F_ROH";next} {print $0, ($3*1000)/aut_len}' ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt > ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.txt
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.txt) > ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.histo 
rclone -v copy ${OUTPUT_DIR}/divStats  --drive-shared-with-me --include "roh_summary_by_RG_L3_Froh.*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/"

awk '{if($5>0.3)print $0}' ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.txt | tr '\t' ',' > ${OUTPUT_DIR}/divStats/roh_high.csv
rclone -v copy ${OUTPUT_DIR}/divStats/roh_high.csv  --drive-shared-with-me "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/"

## Summary stats of all ROH metrics Stratified by the gait type
awk 'BEGIN{FS=OFS="\t";gait["IID"]="gait"}FNR==NR{gait[$2]=$3;next} {if(gait[$1])print $0,gait[$1];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.txt > ${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.txt
INPUT_ROH="${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.txt"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.sumStats.csv"
python scripts/summary_roh.py -i "$INPUT_ROH" -o "$OUTPUT_FILE" -n 4
# Summary saved to ${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.sumStats.csv

## Summary stats of all ROH metrics Stratified by book size for each gait type
awk 'BEGIN{FS=OFS="\t";gait["IID"]="gait"}FNR==NR{gait[$2]=$3;next} {if(gait[$1])print $0,gait[$1];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.txt > ${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt
INPUT_ROH="${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.sumStats.csv"
python scripts/summary_roh.py -i "$INPUT_ROH" -o "$OUTPUT_FILE" -n 4
# Summary saved to ${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.sumStats.csv

## Froh vs ROHshared 
roh_RG="${OUTPUT_DIR}/divStats/roh.L3"
froh="${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt"
for gp in "wholePop" "twoGait" "threeBooksize";do 
    conShare=${roh_RG}.perSample_intersect_${gp}_consensus_${pct}pct.summary.txt  ## the concensus length and % per sample (file for each "rg"), calculated in ROH section 
    output_file="${OUTPUT_DIR}/divStats/Froh_vs_ROHsh_${gp}.png"
    python scripts/roh_plot.py "$conShare" "$froh" "$output_file"
    rclone -v copy "$output_file" --drive-shared-with-me "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/"
    output_prefix="${OUTPUT_DIR}/divStats/normalized_ROHsh_${gp}"
    run_python scripts/roh_histograms.py --metric ratio "$conShare" "$froh" "$output_prefix"
    upload "$output_prefix".histogram.png "Froh/"
    upload "$output_prefix".density.png "Froh/"
    output_prefix2="${OUTPUT_DIR}/divStats/ROHshared_${gp}"
    run_python scripts/roh_histograms.py --metric shared "$conShare" "$froh" "$output_prefix2"
    upload "$output_prefix2".histogram.png "Froh/"
done &> ${OUTPUT_DIR}/divStats/roh_sh.log
## Outputs: (tested in "wholePop" "twoGait" "threeBooksize" BUT the best informative is "twoGait")
## ${OUTPUT_DIR}/divStats/Froh_vs_ROHsh_${gp}.png
## ${OUTPUT_DIR}/divStats/normalized_ROHsh_${gp}
## ${OUTPUT_DIR}/divStats/ROHshared_${gp}

############################################
## x. Nucleotide diversity statistic (pi) -- This section is under development
############################################
## Nucleotide diversity is a population-level metric, the average number of differences between a pair of chromosomes, across all chromosome combinations within the population.
## This is distinct from simply measuring heterozygosity. 
## Variant vs. Invariant Sites: Traditional pi calculations require knowledge of both variant and invariant sites (i.e., sequencing data). 
##    With a SNP array, pi will be overestimated because it ignores the conserved (non-variable) parts of the genome (i.e., it is "SNP-based" diversity rather than a true "genomic" diversity.)



############################################
## 6. Relatedness work
############################################
log "Section 6: Relatedness analysis"
mkdir -p $work_dir/${OUTPUT_DIR}/rep_ROHRM
rohrm_dir="$work_dir/${OUTPUT_DIR}/rep_ROHRM"
group="gait"
phenotypes="${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait"


############################################
## Genomic Relatedness using Plink
############################################
## Create a standard GRM matrix 
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --make-rel square \
    --output-chr 'chrM' --out $rohrm_dir/filtered.LD_prune.GRM_$group

############################################
## ROH-based Relatedness (REPLICATE C++ code of Howard et al.)
############################################
## The C++ code you provided reveals a very specific logic: it relies on phased genotypes (distinguishing between Maternal vs. Paternal alleles), 
## and it calculates relationships based on exact haplotype matching within dynamically defined genomic windows.
# Stage 1: Loads the VCF, converts positions to Megabases (Mb), and splits genotypes into hap1/hap2.
#   To properly implement this, we previously converted the plink files into VCF, select autosomes, bcftools norm to remove duplicate markers, and phase by Beagle.
# Stage 2 (The Geometry Layer): Replicates the logic of the "ROH_Index class" in C++ to identify all valid genomic windows based on physical distance (Mb). 
#   It iterates through every SNP i and finds the window of SNPs that fits within the ROH Cutoff (e.g., 1 Mb). (Sliding Window with a step of 1 SNP)
    # Note: We are not checking for ROH. We are defining the windows only to be used in the next stages.
    # Note2: this implementation is faster than C++ which uses nested loops. Here, we use a single loop with a moving "end" pointer. This reduces the complexity from O(N^2) to O(N).
# Stage 3 (Window Filtering; The Statistical Layer): Discard windows that are too "sparse" (i.e., has too few SNPs) based on the statistical distribution.
    # a. Calculate the Mean and Standard Deviation (SD) of the "Number of SNPs" across all windows.
    # b. Define a threshold: Cutoff = Mean - (User_Threshold * SD).
    # c. Discard any window where NumSNP < Cutoff.
# Stage 4 (The Kernel): Calculate the relationship matrix. This is the hardest part to optimize in Python.
    # a. Extract the vector of genotypes for that window.
    # b. Expand genotypes into haplotypes (Paternal/Maternal).
    # c. Check similarity: If either haplotypes from animal i == either haplotypes from animal j, score 1.0 for any match; otherwise 0.0.
    #.   Note: when $i = j$ (The Diagonal): This function scores the inbreeding coefficient based on ROH: if hap1 == hap2, score 1.0; else 0.0. 
    # d. Accumulate these scores into the global matrix.
    # Note: Again, we optimize the C++ nested loops by using NumPy Broadcasting instead.
# Stage 5 (Normalization): Scale the matrix and save it.
: <<'COMMENT'
Feature,            Standard GRM (VanRaden),                    ROH GRM (Howard et al.)
Input Data,         "Unphased Genotypes (0, 1, 2)",             Phased Haplotypes (0|0, 0|1, 1|0, 1|1)
Resolution,         Single Nucleotide (SNP),                    Multi-Megabase Window (Window)
Matching Logic,     Allele Sharing (IBS),                       Exact String Matching (IBD)
Sensitivity,        Very tolerant of mutation/recombination.,   Very strict. One mismatch breaks the link.
Biological Signal,  Captures Deep/Ancient Relatedness.,         Captures Recent Relatedness.
Diagonal,           Heterozygosity-based Inbreeding.,           ROH-based Inbreeding (FROH​).
COMMENT

############################################
## Compare Genomic Relatedness and ROH-based Relatedness for each window size
## The script outputs "Robust_Matrix_Comparison_Enhanced.png", "Pairwise_Differences.csv", and "Inbreeding_Comparison.csv"
############################################
roh_threshold=$ROH_THRESHOLD_SD
for roh_mb_cutoff in $ROH_CUTOFFS; do
    roh_sd_label="${roh_threshold%.*}SD"
    subfolder="roh_${roh_mb_cutoff%.*}Mb.Threshold_${roh_sd_label}"
    log "Running ROHRM: ${roh_mb_cutoff} Mb cutoff, threshold ${roh_sd_label}"

    run_python $scripts/ROHRM_Creator.py $vcf_filtered.norm.phased.vcf.gz $roh_mb_cutoff $roh_threshold "$rohrm_dir" \
        > $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold.log

    run_python $scripts/analysis_comparison.py \
        $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold \
        $rohrm_dir/filtered.LD_prune.GRM_$group $phenotypes "$rohrm_dir"

    # center column 7 (Difference) around the column's mean and save to new column centered_Kinship_diff
    awk -F, 'BEGIN{OFS=FS=","} NR==1{hdr=$0; next} {sum+=$7; n++; lines[n]=$0; vals[n]=$7} \
             END{ mean = (n?sum/n:0); print hdr, "centered_Kinship_diff"; \
                  for(i=1;i<=n;i++) printf "%s%s%.8f\n", lines[i], OFS, vals[i]-mean }' \
        $rohrm_dir/Pairwise_Differences.csv > $rohrm_dir/Pairwise_Differences.csv.tmp \
        && mv $rohrm_dir/Pairwise_Differences.csv.tmp $rohrm_dir/Pairwise_Differences.csv

    # generate histograms of kinship distributions
    awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                        END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_Std.histo
    awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                        END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_ROH.histo
    awk -v size=0.01 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                        END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_diff.histo

    # extra pair analysis for the primary ROH cutoff only
    if [[ "$roh_mb_cutoff" == "$PRIMARY_ROH_MB" ]]; then
        ## animal pairs with high positive kinship difference (ROH-based kinship > standard kinship)
        awk 'BEGIN{FS=","} NR==1{print;next}{if($8>0.1) print}' \
            $rohrm_dir/Pairwise_Differences.csv > ${OUTPUT_DIR}/divStats/high_positive_kinship_diff.csv
        ## animal pairs with high negative kinship difference (standard kinship > ROH-based kinship)
        awk 'BEGIN{FS=","} NR==1{print;next}{if($8<-0.15) print}' \
            $rohrm_dir/Pairwise_Differences.csv > ${OUTPUT_DIR}/divStats/high_negative_kinship_diff.csv
        ## average centered kinship difference between and within gait groups
        awk 'BEGIN{FS=","} /Trotter/ && /Pacer/ {sum+=$8; n++} END{print "Ave diff Trotter-Pacer:", sum/n}' $rohrm_dir/Pairwise_Differences.csv
        awk 'BEGIN{FS=","} /Trotter/ && !/Pacer/ {sum+=$8; n++} END{print "Ave diff Trotter-Trotter:", sum/n}' $rohrm_dir/Pairwise_Differences.csv
        awk 'BEGIN{FS=","} !/Trotter/ && /Pacer/ {sum+=$8; n++} END{print "Ave diff Pacer-Pacer:", sum/n}' $rohrm_dir/Pairwise_Differences.csv
    fi

    # upload results
    upload $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png "Relatedness/$subfolder/"
    upload $rohrm_dir/Pairwise_Differences.csv "Relatedness/$subfolder/"
    upload $rohrm_dir/Inbreeding_Comparison.csv "Relatedness/$subfolder/"

    # move results to a per-cutoff subfolder
    mkdir -p $rohrm_dir/$subfolder
    mv $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png $rohrm_dir/Inbreeding_Comparison.csv \
       $rohrm_dir/Pairwise_Differences.* $rohrm_dir/$subfolder/
done


## Temp ###########################################
## Compare Inbreeding Coefficients from standard GRM and ROH-based GRM vs Heterozygosity and COIfficient of Inbreeding (COI)
############################################
## compare with het and coi
## Rscript that plots the correlation between KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
RM_diag="${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD/Inbreeding_Comparison.csv" ## to read D_STD (comparable to COI "F" measure in 1) and D_ROH (comparable to F_ROH measured in 5)
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"                ## to read O(HET), E(HET), and F_SNP
Froh_stats="${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.txt"          ## to read F_ROH measured in 5
out_prefix="${OUTPUT_DIR}/divStats/coi_Froh_rmDiag_correlation"
Rscript scripts/correlation_plot.R --mode froh $RM_diag $het_stats $Froh_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me


#roh_RG="${OUTPUT_DIR}/divStats/roh.L3"; rg="wholePop"; pct=25; conShare=${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt;
roh_RG="${OUTPUT_DIR}/divStats/roh.L3"; rg="twoGait"; pct=25; conShare=${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt;
out_prefix2="${OUTPUT_DIR}/divStats/coi_Froh_rmDiag_conShare_correlation"
Rscript scripts/correlation_plot.R --mode froh-cons $RM_diag $het_stats $Froh_stats $conShare $out_prefix2
rclone -v copy $out_prefix2.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## focus on D_STD vs ROH_shared
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$3;next}{print $0,a[$1]}' <(cat $conShare | sed 's/Percent_of_Consensus_ROH/ROH_sh/') <(cat $RM_diag | tr ',' '\t') > ${OUTPUT_DIR}/divStats/rmdiag_conShare

awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize ${OUTPUT_DIR}/divStats/rmdiag_conShare > ${OUTPUT_DIR}/divStats/rmdiag_conShare_wBooksize
input_file="${OUTPUT_DIR}/divStats/rmdiag_conShare_wBooksize"
Rscript $scripts/plot_correlation_withColorsAndShapes.R "$input_file" ROH_sh D_STD Phenotype Book_Size
rclone -v copy ${OUTPUT_DIR}/divStats/correlation_plot_ROH_sh_vs_D_STD_doubleAnn.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## focus on F_ROH vs D_ROH
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$5;next}{print $0,a[$1]}' $Froh_stats <(cat $RM_diag | tr ',' '\t') > ${OUTPUT_DIR}/divStats/rmdiag_Froh
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize ${OUTPUT_DIR}/divStats/rmdiag_Froh > ${OUTPUT_DIR}/divStats/rmdiag_Froh_wBooksize
input_file="${OUTPUT_DIR}/divStats/rmdiag_Froh_wBooksize"
Rscript $scripts/plot_correlation_withColorsAndShapes.R "$input_file" F_ROH D_ROH Phenotype Book_Size
rclone -v copy ${OUTPUT_DIR}/divStats/correlation_plot_F_ROH_vs_D_ROH_doubleAnn.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## focus on F_SNP vs D_ROH
awk 'BEGIN{FS=OFS="\t"}NR==1{a[$2]="F_SNP";next}NR==FNR{a[$2]=$8;next}{print $0,a[$1]}' $het_stats <(cat $RM_diag | tr ',' '\t') > ${OUTPUT_DIR}/divStats/rmdiag_Fsnp
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize ${OUTPUT_DIR}/divStats/rmdiag_Fsnp > ${OUTPUT_DIR}/divStats/rmdiag_Fsnp_wBooksize
input_file="${OUTPUT_DIR}/divStats/rmdiag_Fsnp_wBooksize"
Rscript $scripts/plot_correlation_withColorsAndShapes.R "$input_file" F_SNP D_ROH Phenotype Book_Size
rclone -v copy ${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_D_ROH_doubleAnn.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## focus on F_SNP vs F_ROH
awk 'BEGIN{FS=OFS="\t"}NR==1{a[$2]="F_SNP";next}NR==FNR{a[$2]=$8;next}{print $0,a[$1]}' $het_stats ${OUTPUT_DIR}/divStats/rmdiag_Froh > ${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize ${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp > ${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp_wBooksize
input_file="${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp_wBooksize"
Rscript $scripts/plot_correlation_withColorsAndShapes.R "$input_file" F_SNP F_ROH Phenotype Book_Size
rclone -v copy ${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_F_ROH_doubleAnn.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me


## Resume relatedness work ############################################
############################################
## KING-robust kinship estimator
############################################
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --make-king-table 'counts' 'cols=+ibs1' \
    --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.king_$group

kingkin="${OUTPUT_DIR}/divStats/filtered.LD_prune.king_$group.kin0"
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($10/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 $kingkin) > ${kingkin%.kin0}.histo
rclone -v copy $kingkin "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
rclone -v copy ${kingkin%.kin0}.histo "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## Likely first-degree relations (maybe we neeed to increase this cut-off for such an inbreed population)
head -n 1 $kingkin > ${OUTPUT_DIR}/divStats/related && tail -n +2 $kingkin | sort -grk10,10 | awk '{if($10>0.177)print}' >> ${OUTPUT_DIR}/divStats/related
grep "Trotter" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group | cut -f2 | grep -Fwf - ${OUTPUT_DIR}/divStats/related > ${OUTPUT_DIR}/divStats/related_Trotter
grep "Pacer" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group | cut -f2 | grep -Fwf - ${OUTPUT_DIR}/divStats/related > ${OUTPUT_DIR}/divStats/related_Pacer
############################################
## calc IBS
############################################
## IBS1= HET1_HOM2 + HET2_HOM1 
## IBS2= N_SNPs - (HETHET + IBS0 + IBS1)
## IBS = (2*IBS2 + IBS1) / (2*N_SNPs)
#awk 'BEGIN{FS=OFS="\t"}NR==1{print $0,"IBS";next}{print $0,(2*$7+$8+$9)/(2*$5)}' $kingkin > ${kingkin}.withIBS
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0,"IBS";next}{ibs1=$8+$9;ibs2=$5-($6+$7+ibs1);print $0,(2*ibs2+ibs1)/(2*$5)}' $kingkin > ${kingkin}.withIBS

## useless
## Plot the correlation between KING-robust kinship and IBS
Rscript $scripts/plot_correlation.R ${kingkin}.withIBS KINSHIP IBS
rclone -v copy ${OUTPUT_DIR}/divStats/correlation_plot_KINSHIP_vs_IBS.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
#Plot saved as: correlation_plot_KINSHIP_vs_IBS.png 
#Correlation (Pearson): 0.392 

############################################
## PCA-based pairwise Euclidean distance
############################################
for prefix in ${OUTPUT_DIR}/divStats/filtered.LD_prune{\.,\.Trotter\.,\.Pacer\.}pca ;do 
    awk '
    BEGIN {FS=OFS="\t"}
    $1!="#FID" && NF>3 {
        n++
        fid[n]=$1
        iid[n]=$2
        for(c=3;c<=NF;c++) pc[n,c-2] = $c
        npcs = NF - 2
    }
    END {
        print "FID1","IID1","FID2","IID2","PCA_EUCLIDEAN_DIST","DIST_KINSHIP"
        for(i=1;i<=n;i++)
            for(j=1;j<i;j++) {
                dist2=0
                for(c=1;c<=npcs;c++) {
                    d = pc[i,c] - pc[j,c]
                    dist2 += d*d
                }
                dist = sqrt(dist2)
                printf "%s\t%s\t%s\t%s\t%.8f\t%.8f\n",
                    fid[i], iid[i], fid[j], iid[j], dist, exp(-dist2/2)
            }
    }
    ' "$prefix.eigenvec" > "$prefix.pca_pairwise_euclidean.dist"
done

## Useless
## Merge PCA-based Euclidean distance with KING-robust kinship + IBS
kingkin_wIBS=${kingkin}.withIBS
for prefix in ${OUTPUT_DIR}/divStats/filtered.LD_prune{\.,\.Trotter\.,\.Pacer\.}pca ;do 
    euclDist="$prefix.pca_pairwise_euclidean.dist"
    out_file="$prefix.pca_pairwise_euclidean.dist.withKIN0"
    awk 'BEGIN {FS=OFS="\t"} FNR == NR {
        if ($1 ~ /^#/) next
        # canonical key for pair
        key = ($1":"$2 < $3":"$4) ?
            $1":"$2"|" $3":"$4 :
            $3":"$4"|" $1":"$2
        # store columns 5+ only
        kin0_extra = ""
        for (i = 5; i <= NF; i++)
            kin0_extra = kin0_extra OFS $i

        kin0_data[key] = substr(kin0_extra, 2)   # remove leading OFS
        next
    }
    FNR == 1 {
        print "FID1","IID1","FID2","IID2",
            "PCA_EUCLIDEAN_DIST","KINSHIP_KING_PCA",
            "NSNP","HETHET","IBS0","HET1_HOM2","HET2_HOM1","KINSHIP_PLINK","IBS"
        next
    }
    {
        # build canonical key
        key = ($1":"$2 < $3":"$4) ?
            $1":"$2"|" $3":"$4 :
            $3":"$4"|" $1":"$2

        # fetch .kin0 info (5+ columns)
        extra = (key in kin0_data ? kin0_data[key] : "NA")
        print $1,$2,$3,$4,$5,$6,extra
    }
    ' "$kingkin_wIBS" "$euclDist" > "$out_file"
done

## useless
## Plot the correlation between PCA-based Euclidean distance and KING-robust kinship
for prefix in ${OUTPUT_DIR}/divStats/filtered.LD_prune{\.,\.Trotter\.,\.Pacer\.}pca ;do 
    out_file="$prefix.pca_pairwise_euclidean.dist.withKIN0"
    Rscript $scripts/plot_correlation.R "$out_file" PCA_EUCLIDEAN_DIST KINSHIP_PLINK
    mv ${OUTPUT_DIR}/divStats/correlation_plot_PCA_EUCLIDEAN_DIST_vs_KINSHIP_PLINK.png "$prefix.correlation_plot_PCA_EUCLIDEAN_DIST_vs_KINSHIP_PLINK.png"
    rclone -v copy "$prefix.correlation_plot_PCA_EUCLIDEAN_DIST_vs_KINSHIP_PLINK.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
done

############################################
## Plot the correlation among ROH-based relatedness, KING-robust kinship, and PCA-based Euclidean distance
############################################
RMs="${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD/Pairwise_Differences.csv" ## file1 has Kinship_Std and Kinship_ROH columns
kingkin_wIBS="${OUTPUT_DIR}/divStats/filtered.LD_prune.king_$group.kin0.withIBS" ## file2 has KINSHIP and IBS columns
for pop in "wholePop" "Trotter" "Pacer";do
    pca_prefix=$(echo ${OUTPUT_DIR}/divStats/filtered.LD_prune.$pop.pca | sed 's/wholePop\.//')
    euclDist="$pca_prefix.pca_pairwise_euclidean.dist" ## file3 has PCA_EUCLIDEAN_DIST column
    out_prefix="${OUTPUT_DIR}/divStats/$pop.relatedness_correlation"
    Rscript scripts/correlation_plot.R --mode pairwise $RMs $kingkin_wIBS $euclDist $out_prefix
    rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
done
cat ${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD/Pairwise_Differences.csv | sed 's/ID/IID/g' | awk 'BEGIN{FS=",";OFS="\t";}{a[1]=$1;a[2]=$2;asort(a);print a[1],a[2],$5,$6}' > ${OUTPUT_DIR}/divStats/tmp_kin1
cat ${OUTPUT_DIR}/divStats/filtered.LD_prune.king_$group.kin0.withIBS | awk 'BEGIN{FS=OFS="\t";}{a[1]=$2;a[2]=$4;asort(a);print a[1],a[2],$10,$11}' > ${OUTPUT_DIR}/divStats/tmp_kin2
awk 'BEGIN{FS=OFS="\t";}NR==FNR{a[$1 FS $2]=$0;next}{if(a[$1 FS $2])print a[$1 FS $2],$3,$4}' ${OUTPUT_DIR}/divStats/tmp_kin1 ${OUTPUT_DIR}/divStats/tmp_kin2 > ${OUTPUT_DIR}/divStats/merged_kin
head -n 1 ${OUTPUT_DIR}/divStats/merged_kin > ${OUTPUT_DIR}/divStats/merged_kin_sorted_top && tail -n +2 ${OUTPUT_DIR}/divStats/merged_kin | sort -grk5,5 | awk '{if($5>0.1)print}' >> ${OUTPUT_DIR}/divStats/merged_kin_sorted_top

########################################################
