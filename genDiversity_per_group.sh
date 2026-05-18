#!/usr/bin/env bash
# Per-group metrics stage for the genDiv pipeline.
# Invoked once per $rg ∈ {wholePop, Trotter, Pacer} by genDiversity.sh.
# Produces all per-group reference files consumed by GPA's downstream report
# pipeline, plus the cross-method plots that previously ran per-group inline.
#
# Section order mirrors the original monolithic pipeline:
#   1. afreq (per-group, feeds --read-freq below)
#   2. PCA + overlays
#   3. FST of book-size subgroups within this gait (Trotter / Pacer only)
#   4. Heterozygosity (F_SNP) via --read-freq
#   5. PCA COI overlay (reads the per-group .het from step 4)
#   6. bcftools roh + L1/L2/L3 filter on the group-subset phased VCF
#   7. Per-base consensus ROH (nested over book-size for Trotter / Pacer)
#   8. F_ROH summary (roh_summary_by_RG_L3_Froh.${rg}.txt)
#   9. GRM + ROHRM + analysis_comparison → Inbreeding_Comparison / Pairwise_Differences
set -eo pipefail
source "$(dirname "$0")/genDiversity_common.sh"

rg="${1:?usage: $0 <rg>  (rg must be one of wholePop, Trotter, Pacer)}"
case "$rg" in
    wholePop|Trotter|Pacer) : ;;
    *) echo "ERROR: invalid rg='$rg'; must be wholePop, Trotter, or Pacer" >&2; exit 2 ;;
esac

log "Per-group stage: rg=${rg}"

## Canonical paths produced by genDiversity_shared.sh. Re-declared here because
## subshell variables don't propagate back from shared.sh's invocation.
pl1_filtered="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered"
vcf_filtered="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.vcf"
pl1_pruned="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune"
vcf_pruned="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.LD_prune.vcf"
phenotypes="${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait"
docs="$(pwd)/${OUTPUT_DIR}/Miscellaneous_documents_standardbred"
samples_rg="${OUTPUT_DIR}/preprocess/samples.${rg}.txt"
aut_len=$(cat "${OUTPUT_DIR}/divStats/effective_autosomal_genome_length.txt")

group="gait"

## Strict GPA naming: every per-group file has ".${rg}." inserted before
## its extension, including for wholePop. Consumer update contract in
## MIGRATION.md.
rg_tag=".${rg}"

mkdir -p "${OUTPUT_DIR}/divStats" "${OUTPUT_DIR}/LD_pruned"

##########################################
## 1. Per-group allele frequency on the LD-pruned SNP set
##########################################
pruned_afreq="${OUTPUT_DIR}/LD_pruned/pruned.${rg}.afreq"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$samples_rg" \
    --freq \
    --output-chr 'chrM' --out "${OUTPUT_DIR}/LD_pruned/pruned.${rg}"
rclone -v copy "$pruned_afreq" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Ae/" --drive-shared-with-me

##########################################
## 2. PCA + overlays
##########################################
## PCAs are called "loadings" because they represent the weights or coefficients that 
## determine how much each original variable "loads" onto or contributes to a specific PC.
pca_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune${rg_tag}.pca"
if [[ "$rg" == "wholePop" ]]; then
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
           --real-ref-alleles --autosome --pca 'allele-wts' \
           --read-freq "$pruned_afreq" \
           --output-chr 'chrM' --out "$pca_prefix"
    n_pcs=6
else
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
           --keep "$samples_rg" --autosome --pca \
           --read-freq "$pruned_afreq" \
           --output-chr 'chrM' --out "$pca_prefix"
    n_pcs=3
fi

rclone -v copy "${OUTPUT_DIR}/divStats" --include "$(basename "$pca_prefix").eigen*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

Rscript -e 'args=(commandArgs(TRUE));'\
'val <- read.table(paste(args[1],"eigenval",sep="."));'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = args[2]);'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "principal Component", ylab = "Variance explained in %");'\
'dev.off();' "$pca_prefix" "${OUTPUT_DIR}/divStats/Var_PCs${rg_tag}.jpg"
rclone -v copy "${OUTPUT_DIR}/divStats/Var_PCs${rg_tag}.jpg" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# Color by Book Size (all groups)
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize" "$pca_prefix.eigenvec" > "$pca_prefix.eigenvec.wBook_Size"
eigenvec_suffix="wBook_Size"; color_column="Book_Size"
out_png="${OUTPUT_DIR}/divStats/pca_plot_BookSize.${rg}.png"
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" "$n_pcs" factor
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

# wholePop-only: Sex overlay, Gait overlay, PC outlier extractions
if [[ "$rg" == "wholePop" ]]; then
    # Sex overlay
    awk 'BEGIN{FS=OFS="\t";a["IID"]="sex"}NR==FNR{if($5==1)a[$2]="male";else a[$2]="female";next}{print $0,a[$2]}' "$pl1_pruned.fam" "$pca_prefix.eigenvec" > "$pca_prefix.eigenvec.wSex"
    eigenvec_suffix="wSex"; color_column="sex"; out_png="${OUTPUT_DIR}/divStats/pca_plot_sex.png"
    Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" 6 factor
    rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

    # Gait overlay
    awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait" "$pca_prefix.eigenvec" > "$pca_prefix.eigenvec.wGait"
    eigenvec_suffix="wGait"; color_column="Gait"; out_png="${OUTPUT_DIR}/divStats/pca_plot_Gait.png"
    Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" 6 factor
    rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

    ## Identify Trotter samples segregating on PC2 (actually PC4 after highly-related exclusion)
    cat "$pca_prefix.eigenvec" | awk 'BEGIN{FS=OFS="\t"}{if($6>0.1)print $2}' | grep -Fwf - "$docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv" > "${OUTPUT_DIR}/divStats/Trotters_segregating_on_PC2.csv" || true
    rclone -v copy "${OUTPUT_DIR}/divStats/Trotters_segregating_on_PC2.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
    ## Pacers co-segregating with Trotters on PC1
    cat "$pca_prefix.eigenvec" | awk 'BEGIN{FS=OFS="\t"}{if($3<0)print $2}' | grep -Fwf - "$docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv" | grep "Pacer" > "${OUTPUT_DIR}/divStats/Pacers_cosegregating_withTrotters_on_PC1.csv" || true
    rclone -v copy "${OUTPUT_DIR}/divStats/Pacers_cosegregating_withTrotters_on_PC1.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
    ## Trotters co-segregating with Pacers on PC1
    cat "$pca_prefix.eigenvec" | awk 'BEGIN{FS=OFS="\t"}{if($3>0)print $2}' | grep -Fwf - "$docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv" | grep "Trotter" > "${OUTPUT_DIR}/divStats/Trotters_cosegregating_withPacers_on_PC1.csv" || true
    rclone -v copy "${OUTPUT_DIR}/divStats/Trotters_cosegregating_withPacers_on_PC1.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
fi

########################################################
## 1. Effective number of alleles (\(A_{e}\))
## A_e represents the number of equally frequent alleles required to achieve the 
## same level of expected heterozygosity (\(H_{e}\)) observed in a population

## Implementation note: The A_e is still implemented in the old design using wholePop to calculate all gait subpopulations stats
########################################################
## Formula
## A_e = 1/Σ p_i^2 per SNP
## where (p_i) is the frequency of the (i^{th}) allele

## Example Calculation
## Suppose a single locus has three alleles with the observed frequencies (0.6, 0.3, 0.1) in a population.
## 1. Calculate the squared frequencies: 0.36, 0.09, and 0.01
## 2. Calculate \(A_{e}\): 1/(0.36 + 0.09 + 0.01) = 1/0.46 = 2.17
## This result means that although there are 3 distinct alleles, the population's genetic diversity is equivalent to a population with only 2.17 equally frequent alleles.

## Calculate \(A_{e}\) for each SNP
if [[ "$rg" == "wholePop" ]]; then
    ## Calculate \(A_{e}\) for each SNP
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --freq \
        --out "$pl1_pruned.freq_stats"
    awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.afreq" > "$pl1_pruned.freq_stats.wholePop.afreq.Ae"
    awk -v pop="wholePop" 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
        { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"pop, mean_Ae, "SD_Ae_in_"pop, sd_Ae } }' "$pl1_pruned.freq_stats.wholePop.afreq.Ae"
    #Mean_Ae_in_wholePop     1.5616  SD_Ae_in_wholePop       0.321994

    ## Calculate \(A_{e}\) for each SNP per gait subpopulation
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
fi

##########################################
## 2. Fst between subpopulations (genders, gait types, and book sizes)
##########################################
## The fixation index can range from 0 to 1, where 0 means complete sharing of genetic material and 1 means no sharing.
## For values equal to 1(meaning no sharing), scientists say that the populations are fixed.
## Effects of marker type and filtering criteria on QST-FST comparisons: https://pmc.ncbi.nlm.nih.gov/articles/PMC6894560/
if [[ "$rg" == "wholePop" ]]; then
    for fst_group in sex gait bookSize; do
        plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
            --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$fst_group \
            --fst 'PHENO1' 'blocksize=2000' \
            --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_$fst_group
    done

    find ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_*.summary -maxdepth 1 -type f | grep -v "\.x\." | xargs cat > ${OUTPUT_DIR}/divStats/autosomal.fst.summary
    rclone -v copy ${OUTPUT_DIR}/divStats/autosomal.fst.summary "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me
fi

##########################################
## 3. FST of book-size subgroups within this gait (Trotter / Pacer only)
##########################################
if [[ "$rg" != "wholePop" ]]; then
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --keep "$samples_rg" \
        --pheno "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize" \
        --fst 'PHENO1' 'blocksize=2000' \
        --output-chr 'chrM' --out "${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_bookSize.${rg}"
    rclone -v copy "${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_bookSize.${rg}.fst.summary" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me
fi

##########################################
## 4. Expected and observed heterozygosity and inbreeding coefficient (COI) -- aka (F_SNP) 
## using group-specific AF via --read-freq ($pruned_afreq is now $rg dependanat)
##########################################
het_rg_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.${rg}"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$samples_rg" \
    --het 'cols=fid,hom,het,nobs,f' \
    --read-freq "$pruned_afreq" \
    --output-chr 'chrM' --out "$het_rg_prefix"
rclone -v copy "${het_rg_prefix}.het" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/het_and_COI/" --drive-shared-with-me

##########################################
## 4b. Whole-pop .het × gait / gait_bookSize stratified summaries (wholePop only)
##########################################
## Joins the whole-pop per-sample .het (just produced above) with the global
## gait / gait_bookSize metadata from shared.sh, then runs summary_het.py on
## each stratification. Only uses wholePop inputs, so runs once.
if [[ "$rg" == "wholePop" ]]; then
    awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' \
        "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait" "${het_rg_prefix}.het" \
        > "${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait"
    python scripts/summary_het.py \
        -i "${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait" \
        -o "${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait.sumStats.csv"

    awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' \
        "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize" "${het_rg_prefix}.het" \
        > "${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize"
    python scripts/summary_het.py \
        -i "${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize" \
        -o "${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize.sumStats.csv"

    rclone -v copy "${OUTPUT_DIR}/divStats" --include "filtered.LD_prune.het_stats.het.wGait*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/het_and_COI/" --drive-shared-with-me
fi

##########################################
## 5. PCA COI overlay (uses this group's own .het)
##########################################
awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 "${het_rg_prefix}.het") "$pca_prefix.eigenvec" > "$pca_prefix.eigenvec.wCOI"
eigenvec_suffix="wCOI"; color_column="COI"
out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.${rg}.png"
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" "$n_pcs" numeric
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

##########################################
## 6. ROH using bcftools/roh (Filtered dataset without LD pruning) + L1/L2/L3 filter 
## Implementation note: Calc using group-subset VCF so AF is group-specific
##########################################
## For wholePop, we use the whole-pop phased VCF directly; for Trotter / Pacer
## We create a group-subset phased VCF that is ALSO reused by ROHRM in §9.
if [[ "$rg" == "wholePop" ]]; then
    group_vcf="${vcf_filtered}.norm.phased.vcf.gz"
else
    group_vcf="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.${rg}.vcf.gz"
    bcftools view -S <(cut -f2 "$samples_rg") --force-samples \
        "${vcf_filtered}.norm.phased.vcf.gz" -Oz -o "$group_vcf"
    bcftools index -t "$group_vcf"
fi

## Per-group tabix-indexed AF table consumed by downstream GPA's `bcftools roh
## --AF-file` on 1-2 animal mate-pair VCFs (which can't --estimate-AF
## themselves). Built on the same $group_vcf used for the upstream ROH call
## below, so the AF basis matches what the reference-distribution ROH run sees.
## For wholePop this is numerically equivalent to GPA's legacy self-built
## freqs.tab.gz; for Trotter/Pacer it is new per-group content.
freqs_prefix="${OUTPUT_DIR}/divStats/freqs.${rg}"
bcftools +fill-tags "$group_vcf" -- -t AF \
    | bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
    | bgzip -c > "${freqs_prefix}.tab.gz"
tabix -s1 -b2 -e2 "${freqs_prefix}.tab.gz"
rclone -v copy "${freqs_prefix}.tab.gz"     "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
rclone -v copy "${freqs_prefix}.tab.gz.tbi" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

## ROH calls driven by the group's own AF (via --estimate-AF on the subset VCF).
bcftools roh -G30 --estimate-AF - "$group_vcf" -o "${OUTPUT_DIR}/divStats/roh_out.${rg}.txt"
grep -E "^# RG|^RG" "${OUTPUT_DIR}/divStats/roh_out.${rg}.txt" > "${OUTPUT_DIR}/divStats/roh_out_RG.${rg}.txt"
awk '/^#/ || $6 >= 1000000' "${OUTPUT_DIR}/divStats/roh_out_RG.${rg}.txt" > "${OUTPUT_DIR}/divStats/roh.L1.${rg}.txt"
awk '/^#/ || $7 >= 50'      "${OUTPUT_DIR}/divStats/roh.L1.${rg}.txt"    > "${OUTPUT_DIR}/divStats/roh.L2.${rg}.txt"
awk '/^#/ || $8 >= 20'      "${OUTPUT_DIR}/divStats/roh.L2.${rg}.txt"    > "${OUTPUT_DIR}/divStats/roh.L3.${rg}.txt"

## Per-sample L3 counts (pre-F_ROH; F_ROH column appended in §8).
awk 'BEGIN{print "IID\tNSEG\tKB\tKBAVG"} $1=="RG"{n[$2]++; sum[$2]+=$6} END{for (s in n) printf "%s\t%d\t%.2f\t%.2f\n", s, n[s], sum[s]/1000, (sum[s]/1000)/n[s]}' \
    "${OUTPUT_DIR}/divStats/roh.L3.${rg}.txt" > "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.${rg}.txt"

##########################################
## 6b. Whole-pop ROH-vs-het sanity-check correlation plot (wholePop only)
##########################################
## Basic pairplot correlating NSEG / KB / KBAVG (from roh_summary_by_RG_L3,
## just produced above) with O(HET) / E(HET) / F (from §4's .het). Uses only
## wholePop inputs.
if [[ "$rg" == "wholePop" ]]; then
    out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_summary_by_RG_L3"
    Rscript scripts/correlation_plot.R --mode basic \
        "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.${rg}.txt" \
        "${het_rg_prefix}.het" \
        "$out_prefix"
    rclone -v copy "$out_prefix.pairplot.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
fi

##########################################
## 7. Per-base consensus ROH (this group + its book-size subgroups)
##########################################
## A per-base consensus ROH where ≥25% of "group-specific" samples are in ROH filtered by minimum size 500 kb and stratified by gait type
## With and without applying a smoothing function to the per-base coverage data to reduce noise before identifying consensus ROH regions
##
## Subgroups:
##   wholePop → just wholePop
##   Trotter  → Trotter, Trotter_LOW, Trotter_MEDIUM, Trotter_HIGH
##   Pacer    → Pacer,   Pacer_LOW,   Pacer_MEDIUM,   Pacer_HIGH
##
## The book-size subgroups reuse the gait's ROH calls (Trotter's AF for
## Trotter_*, Pacer's AF for Pacer_*) and differ only in sample membership.
roh_RG="${OUTPUT_DIR}/divStats/roh.L3"
# 7.1 Convert RG output → BED
awk 'BEGIN{OFS="\t"} $1=="RG" {print $3, $4-1, $5, $2}' "${roh_RG}.${rg}.txt" > "${roh_RG}.${rg}.bed"
# 7.2 Merge per sample so overlapping calls in one sample aren't double-counted
cut -f4 "${roh_RG}.${rg}.bed" | sort -u | while read S; do
    awk -v s="$S" '$4==s' "${roh_RG}.${rg}.bed" | sort -k1,1 -k2,2n | bedtools merge -i - | awk -v s="$S" 'BEGIN{OFS="\t"}{print $1,$2,$3,s}'
done | sort -k1,1 -k2,2n > "${roh_RG}.merged_per_sample.${rg}.bed"

# Determine the list of subs to process for this rg
if [[ "$rg" == "wholePop" ]]; then
    subs=("wholePop")
else
    subs=("$rg" "${rg}_LOW" "${rg}_MEDIUM" "${rg}_HIGH")
    # Book-size subsets reuse this gait's per-sample bed, filtered to book-size membership.
    # Tolerate grep's exit-1-when-no-match (small/absent subgroups) by truncating the output.
    for sub in "${rg}_LOW" "${rg}_MEDIUM" "${rg}_HIGH"; do
        (grep "$sub" "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize" | cut -f2 \
            | grep -f - "${roh_RG}.merged_per_sample.${rg}.bed" \
            > "${roh_RG}.merged_per_sample.${sub}.bed") || : > "${roh_RG}.merged_per_sample.${sub}.bed"
    done
fi

# 7.3 Per-base coverage + histograms + consensus threshold + smoothing + per-sample intersect
# autosomes.genome (chromosome sizes table) produced once in shared.sh; reuse.
{
for sub in "${subs[@]}"; do
    bed_perSample="${roh_RG}.merged_per_sample.${sub}.bed"
    bedtools genomecov -i "$bed_perSample" -g "${OUTPUT_DIR}/divStats/autosomes.genome" -bg > "${roh_RG}.per_base_coverage.${sub}.bed"
    awk -v size=5 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($4/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                          END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' "${roh_RG}.per_base_coverage.${sub}.bed" > "${roh_RG}.per_base_coverage.${sub}.histo"

    # upload bed files
    # pause for now to save space
    #rclone -v copy ${roh_RG}.per_base_coverage.${sub}.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/freq/" --drive-shared-with-me

    # upload histo files
    rclone -v copy "${roh_RG}.per_base_coverage.${sub}.histo" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/freq/" --drive-shared-with-me

    num_samples=$(cut -f4 "$bed_perSample" | sort -u | wc -l)
    threshold=$(echo "$pct * $num_samples / 100" | bc -l)

    # Pre-smoothing consensus
    awk -v threshold=$threshold 'BEGIN{OFS="\t"} $4 >= threshold {print}' "${roh_RG}.per_base_coverage.${sub}.bed" > "${roh_RG}.consensus_${pct}pct.${sub}.bed"
    bedtools merge -i "${roh_RG}.consensus_${pct}pct.${sub}.bed" -c 4 -o mean | awk -v min_mb=$CONSENSUS_MIN_MB 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=min_mb)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${sub}.bed"

    echo "==== consensus ROH in ≥${pct}% of ${sub} samples BEFORE SMOOTHING ======"
    awk -v rg="$sub" -v nsam="$num_samples" 'BEGIN{OFS=",";maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
      {print rg,"\nNo. of segments","Total length (KB)","Ave. length (KB)","Max % of samples in consensus","Average % of samples in consensus",\
      "\n"NR,sumLen,sumLen/NR,(maxConsen/nsam)*100"%",((sumSamples/NR)/nsam)*100"%"}' "${roh_RG}.consensus_${pct}pct.merged.${sub}.bed"

    # Smooth per-base coverage: average a middle interval if both flanks are ≥ threshold and middle is <.
    awk -v thr="$threshold" 'BEGIN{OFS="\t"} {chr[NR]=$1; st[NR]=$2; en[NR]=$3; cov[NR]=$4} END{
        for(i=1;i<=NR;i++){
            newcov=cov[i]
            if(i>1 && i<NR){
                if(en[i-1]==st[i] && en[i]==st[i+1]){
                    if(cov[i-1] >= thr && cov[i+1] >= thr && cov[i] < thr){
                        newcov = (cov[i-1] + cov[i] + cov[i+1]) / 3
                    }
                }
            }
            printf "%s\t%d\t%d\t%.6f\n", chr[i], st[i], en[i], newcov
        }
    }' "${roh_RG}.per_base_coverage.${sub}.bed" | awk 'BEGIN{OFS="\t"}{$4=$4+0;print}' > "${roh_RG}.per_base_coverage.${sub}.smoothed.bed"
    awk -v threshold="$threshold" 'BEGIN{OFS="\t"} $4 >= threshold {print}' "${roh_RG}.per_base_coverage.${sub}.smoothed.bed" > "${roh_RG}.consensus_${pct}pct.${sub}.smoothed.bed"
    bedtools merge -i "${roh_RG}.consensus_${pct}pct.${sub}.smoothed.bed" -c 4 -o mean | awk -v min_mb=$CONSENSUS_MIN_MB 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=min_mb)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${sub}.smoothed.bed"

    # upload bed files
    # pause for now to save space
    #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${sub}.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
    #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${sub}.smoothed.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

    echo "==== consensus smoothed ROH in ≥${pct}% of ${sub} samples AFTER SMOOTHING ======"
    awk -v rg="$sub" -v nsam="$num_samples" 'BEGIN{OFS=",";maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
     {print rg,"\nNo. of segments","Total length (KB)","Ave. length (KB)","Max % of samples in consensus","Average % of samples in consensus",\
      "\n"NR,sumLen,sumLen/NR,(maxConsen/nsam)*100"%",((sumSamples/NR)/nsam)*100"%"}' "${roh_RG}.consensus_${pct}pct.merged.${sub}.smoothed.bed"
done
} > "${roh_RG}.consensus_${pct}pct.${rg}.summary.txt"
rclone -v copy "${roh_RG}.consensus_${pct}pct.${rg}.summary.txt" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

# Per-sample intersect against each sub's (smoothed) consensus
for sub in "${subs[@]}"; do
    consensus_bed="${roh_RG}.consensus_${pct}pct.merged.${sub}.smoothed.bed"
    consensus_size=$(awk 'BEGIN{sum=0} {sum+=($3-$2)} END {print sum}' "$consensus_bed")
    bed_perSample="${roh_RG}.merged_per_sample.${sub}.bed"
    echo -e "IID\tTotal_ROH_in_Consensus_region(bp)\tPercent_of_Consensus_ROH" > "${roh_RG}.perSample_intersect_${sub}_consensus_${pct}pct.summary.txt"
    cut -f4 "$bed_perSample" | sort -u | while read S; do
        awk -v s="$S" '$4==s' "$bed_perSample" | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b "$consensus_bed" | awk -v s="$S" -v cs="$consensus_size" 'BEGIN{OFS="\t"}{size+=($3-$2)} END {print s, size, (size/cs)*100}'
    done >> "${roh_RG}.perSample_intersect_${sub}_consensus_${pct}pct.summary.txt"
    #rclone -v copy ${roh_RG}.perSample_intersect_${sub}_consensus_${pct}pct.summary.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
done

##########################################
## 7.5 ROH_common: continuous population-autozygosity landscape + per-individual scores
##########################################
## Builds a fixed 100 kb window tiling for the main group, counts distinct
## IIDs per window (n_w), and writes f_w = n_w/N as the "landscape". Then
## scores each individual against that landscape with an exact leave-one-out
## correction: ROH_common,i = mean over w in W_i of (n_w - 1)/(N - 1).
## Length-stratified landscapes (1-3, 3-5, 5-10, >10 Mb) are emitted in parallel.

roh_common_dir="${OUTPUT_DIR}/divStats/roh_common"
mkdir -p "$roh_common_dir"

run_python scripts/roh_common_landscape.py \
    --segments-bed "${roh_RG}.${rg}.bed" \
    --autosomes    "${OUTPUT_DIR}/divStats/autosomes.genome" \
    --window-kb    "$ROH_COMMON_WINDOW_KB" \
    --length-bins  "$ROH_COMMON_LENGTH_BINS" \
    --out-dir      "$roh_common_dir" \
    --rg           "$rg"

run_python scripts/roh_common_individual.py \
    --landscape-dir "$roh_common_dir" \
    --segments-bed  "${roh_RG}.${rg}.bed" \
    --rg            "$rg" \
    --out           "${roh_common_dir}/roh_common.${rg}.tsv"

upload "${roh_common_dir}/roh_common.${rg}.tsv" "ROH/roh_common/"
for variant in "" ".1to3" ".3to5" ".5to10" ".more10"; do
    f="${roh_common_dir}/landscape.${rg}${variant}.tsv"
    [[ -f "$f" ]] && upload "$f" "ROH/roh_common/"
done

##########################################
## 8. F_ROH statistic and F_ROH summary
##########################################
## F_ROH is an inbreeding coefficient based on runs of homozygosity
## Standard practice is to calculate F_{ROH} statistics on all valid ROHs (>1Mb), while restricting "Islands" (signatures of selection) to only the most robust regions.
## per-sample F_ROH = (sum length of ROH for that individual) / (total autosomal genome length).
## "total autosomal genome length" was calculated in shared.sh

awk -v aut_len="$aut_len" 'BEGIN{FS=OFS="\t";}NR==1{print $0,"F_ROH";next} {print $0, ($3*1000)/aut_len}' \
    "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.${rg}.txt" > "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt"
rclone -v copy "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me

##########################################
## 8b. Whole-pop F_ROH histogram, high-F_ROH shortlist, and gait / gait_bookSize
##     stratified F_ROH + NSEG-length-bin summaries (wholePop only)
##########################################
## Only uses wholePop F_ROH + global metadata, so runs once.
if [[ "$rg" == "wholePop" ]]; then
    froh_wholePop="${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt"

    awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 "$froh_wholePop") > "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.histo"
    rclone -v copy "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.histo" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me

    awk '{if($5>0.3)print $0}' "$froh_wholePop" | tr '\t' ',' > "${OUTPUT_DIR}/divStats/roh_high.csv"
    rclone -v copy "${OUTPUT_DIR}/divStats/roh_high.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me

    awk 'BEGIN{FS=OFS="\t";gait["IID"]="gait"}FNR==NR{gait[$2]=$3;next} {if(gait[$1])print $0,gait[$1];else print $0,"undefined";}' \
        "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait" "$froh_wholePop" \
        > "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.txt"
    python scripts/summary_roh.py \
        -i "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.txt" \
        -o "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.sumStats.csv" -n 4
    rclone -v copy "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.sumStats.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me
    python scripts/summary_nseg_bins.py \
        -i "${OUTPUT_DIR}/divStats/roh.L3.${rg}.txt" \
        -f "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.txt" \
        -o "${OUTPUT_DIR}/divStats/roh.L3_NSEGbins_gait.sumStats.csv"
    rclone -v copy "${OUTPUT_DIR}/divStats/roh.L3_NSEGbins_gait.sumStats.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me

    awk 'BEGIN{FS=OFS="\t";gait["IID"]="gait"}FNR==NR{gait[$2]=$3;next} {if(gait[$1])print $0,gait[$1];else print $0,"undefined";}' \
        "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize" "$froh_wholePop" \
        > "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt"
    python scripts/summary_roh.py \
        -i "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt" \
        -o "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.sumStats.csv" -n 4
    rclone -v copy "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.sumStats.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me
    python scripts/summary_nseg_bins.py \
        -i "${OUTPUT_DIR}/divStats/roh.L3.${rg}.txt" \
        -f "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt" \
        -o "${OUTPUT_DIR}/divStats/roh.L3_NSEGbins_gait_bookSize.sumStats.csv"
    rclone -v copy "${OUTPUT_DIR}/divStats/roh.L3_NSEGbins_gait_bookSize.sumStats.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me
fi

##########################################
## 9. Relatedness Matrices (GRM + ROHRM) → G_SNP, D_SNP, G_ROH, D_ROH 
##########################################
## Per-group standard GRM (vanRaden) via plink2 --make-rel + per-group ROHRM (Howard re-implementation). 
## analysis_comparison.py writes Inbreeding_Comparison.csv (IID, D_STD, D_ROH, Phenotype) and
## Pairwise_Differences.csv (ID1, ID2, Pheno1, Pheno2, Kinship_Std, Kinship_ROH, Difference). 

grp_workdir="${OUTPUT_DIR}/rep_ROHRM/perGroup_${rg}"
mkdir -p "$grp_workdir"

## 9a. Per-group standard GRM (computed once, reused across all ROHRM cutoffs)
grm_prefix="${grp_workdir}/filtered.LD_prune.GRM.${rg}"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$samples_rg" \
    --make-rel 'square' \
    --read-freq "$pruned_afreq" \
    --output-chr 'chrM' --out "$grm_prefix"

## 9b. Per-group ROHRM + analysis_comparison at every cutoff in $ROH_CUTOFFS
## Per-cutoff outputs land at rep_ROHRM/roh_${cutoff%.*}Mb.Threshold_${sd}SD/ with a .${rg}. suffix.
## canonical_dir below is the primary-cutoff subfolder, used later in §14
## (cross-method correlation) and §15 (froh-cons correlation plots) — which
## consume only the primary cutoff's Pairwise_Differences / Inbreeding_Comparison.
rohrm_sd="$ROH_THRESHOLD_SD"
roh_sd_label="${rohrm_sd%.*}SD"
canonical_dir="${OUTPUT_DIR}/rep_ROHRM/roh_${PRIMARY_ROH_MB%.*}Mb.Threshold_${roh_sd_label}"

for rohrm_mb in $ROH_CUTOFFS; do
    subfolder="${OUTPUT_DIR}/rep_ROHRM/roh_${rohrm_mb%.*}Mb.Threshold_${roh_sd_label}"
    mkdir -p "$subfolder"
    log "Running ROHRM (rg=${rg}): ${rohrm_mb} Mb cutoff, threshold ${roh_sd_label}"

    rohrm_prefix="${grp_workdir}/ROHRM.rohMinSize_${rohrm_mb}.rohThreshold_${rohrm_sd}"
    run_python "$scripts/ROHRM_Creator.py" "$group_vcf" "$rohrm_mb" "$rohrm_sd" "$grp_workdir" \
        > "${rohrm_prefix}.log"

    run_python "$scripts/analysis_comparison.py" \
        "$rohrm_prefix" \
        "$grm_prefix" \
        "$phenotypes" "$grp_workdir"

    # Center the Difference column (col 7) on its mean, append as centered_Kinship_diff
    awk -F, 'BEGIN{OFS=FS=","} NR==1{hdr=$0; next} {sum+=$7; n++; lines[n]=$0; vals[n]=$7} \
             END{ mean = (n?sum/n:0); print hdr, "centered_Kinship_diff"; \
                  for(i=1;i<=n;i++) printf "%s%s%.8f\n", lines[i], OFS, vals[i]-mean }' \
        "${grp_workdir}/Pairwise_Differences.csv" > "${grp_workdir}/Pairwise_Differences.csv.tmp" \
        && mv "${grp_workdir}/Pairwise_Differences.csv.tmp" "${grp_workdir}/Pairwise_Differences.csv"

    # Kinship distribution histograms
    awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                          END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 "${grp_workdir}/Pairwise_Differences.csv") > "${grp_workdir}/Pairwise_Differences.Kinship_Std.histo"
    awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                          END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 "${grp_workdir}/Pairwise_Differences.csv") > "${grp_workdir}/Pairwise_Differences.Kinship_ROH.histo"
    awk -v size=0.01 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                          END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 "${grp_workdir}/Pairwise_Differences.csv") > "${grp_workdir}/Pairwise_Differences.Kinship_diff.histo"

    # Primary cutoff only: per-group high-positive / high-negative kinship pair shortlists + gait-mix averages
    if [[ "$rohrm_mb" == "$PRIMARY_ROH_MB" ]]; then
        awk 'BEGIN{FS=","} NR==1{print;next}{if($8>0.1) print}' \
            "${grp_workdir}/Pairwise_Differences.csv" > "${OUTPUT_DIR}/divStats/high_positive_kinship_diff.${rg}.csv"
        awk 'BEGIN{FS=","} NR==1{print;next}{if($8<-0.15) print}' \
            "${grp_workdir}/Pairwise_Differences.csv" > "${OUTPUT_DIR}/divStats/high_negative_kinship_diff.${rg}.csv"
        awk 'BEGIN{FS=","} /Trotter/ && /Pacer/ {sum+=$8; n++} END{if(n)print "Ave diff Trotter-Pacer:", sum/n}' "${grp_workdir}/Pairwise_Differences.csv"
        awk 'BEGIN{FS=","} /Trotter/ && !/Pacer/ {sum+=$8; n++} END{if(n)print "Ave diff Trotter-Trotter:", sum/n}' "${grp_workdir}/Pairwise_Differences.csv"
        awk 'BEGIN{FS=","} !/Trotter/ && /Pacer/ {sum+=$8; n++} END{if(n)print "Ave diff Pacer-Pacer:", sum/n}' "${grp_workdir}/Pairwise_Differences.csv"
    fi

    # Move csvs to the per-cutoff canonical dir with .${rg}. suffix; rename
    # remaining artifacts (matrix comparison png + histos) to carry .${rg}.
    mv "${grp_workdir}/Inbreeding_Comparison.csv" "${subfolder}/Inbreeding_Comparison.${rg}.csv"
    mv "${grp_workdir}/Pairwise_Differences.csv"  "${subfolder}/Pairwise_Differences.${rg}.csv"
    mv "${grp_workdir}/Pairwise_Differences.Kinship_Std.histo"  "${subfolder}/Pairwise_Differences.Kinship_Std.${rg}.histo"
    mv "${grp_workdir}/Pairwise_Differences.Kinship_ROH.histo"  "${subfolder}/Pairwise_Differences.Kinship_ROH.${rg}.histo"
    mv "${grp_workdir}/Pairwise_Differences.Kinship_diff.histo" "${subfolder}/Pairwise_Differences.Kinship_diff.${rg}.histo"
    if [[ -f "${grp_workdir}/Robust_Matrix_Comparison_Enhanced.png" ]]; then
        mv "${grp_workdir}/Robust_Matrix_Comparison_Enhanced.png" "${subfolder}/Robust_Matrix_Comparison_Enhanced.${rg}.png"
    fi

    # Upload per-cutoff outputs under Relatedness/<subfolder-tail>/ with .${rg}. in filenames
    subfolder_tail="$(basename "$subfolder")"
    if [[ -f "${subfolder}/Robust_Matrix_Comparison_Enhanced.${rg}.png" ]]; then
        upload "${subfolder}/Robust_Matrix_Comparison_Enhanced.${rg}.png" "Relatedness/${subfolder_tail}/"
    fi
    upload "${subfolder}/Inbreeding_Comparison.${rg}.csv" "Relatedness/${subfolder_tail}/"
    upload "${subfolder}/Pairwise_Differences.${rg}.csv"  "Relatedness/${subfolder_tail}/"
done

##########################################
## 10. Per-group KING-robust kinship + IBS + related
##########################################
## Computes KING on this group's own samples (--keep samples.${rg}.txt). Each
## group gets its own pair-kinship space; for wholePop, --keep is the full
## cohort so the output is numerically equivalent to the pre-refactor whole-pop
## KING table. Produces: .${rg}.kin0 + .histo, related.${rg} (first-degree
## relatives at KING>0.177), and .${rg}.kin0.withIBS (IBS-augmented kin0).
king_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.king_${group}.${rg}"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$samples_rg" \
    --make-king-table 'counts' 'cols=+ibs1' \
    --output-chr 'chrM' --out "$king_prefix"

kingkin="${king_prefix}.kin0"
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($10/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' <(tail -n+2 "$kingkin") > "${kingkin%.kin0}.histo"
rclone -v copy "$kingkin" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
rclone -v copy "${kingkin%.kin0}.histo" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## Likely first-degree relations (per-group)
head -n 1 "$kingkin" > "${OUTPUT_DIR}/divStats/related.${rg}"
tail -n +2 "$kingkin" | sort -grk10,10 | awk '{if($10>0.177)print}' >> "${OUTPUT_DIR}/divStats/related.${rg}"

## IBS augmentation of kin0: IBS1 = HET1_HOM2 + HET2_HOM1; IBS2 = N_SNPs - (HETHET + IBS0 + IBS1); IBS = (2*IBS2 + IBS1) / (2*N_SNPs).
kingkin_wIBS="${kingkin}.withIBS"
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0,"IBS";next}{ibs1=$8+$9;ibs2=$5-($6+$7+ibs1);print $0,(2*ibs2+ibs1)/(2*$5)}' "$kingkin" > "$kingkin_wIBS"

##########################################
## 11. PCA-based pairwise Euclidean distance
##########################################
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
' "$pca_prefix.eigenvec" > "$pca_prefix.pca_pairwise_euclidean.dist"

##########################################
## 14. Cross-method correlation plot (ROHRM vs KING vs PCA, per group)
##########################################
euclDist="$pca_prefix.pca_pairwise_euclidean.dist"
out_prefix="${OUTPUT_DIR}/divStats/${rg}.relatedness_correlation"
Rscript scripts/correlation_plot.R --mode pairwise \
    "${canonical_dir}/Pairwise_Differences.${rg}.csv" \
    "$kingkin_wIBS" \
    "$euclDist" \
    "$out_prefix"
rclone -v copy "$out_prefix.pairplot.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

##########################################
## 15. F_SNP / F_ROH / D_STD / D_ROH / ROH_sh correlation plots
##########################################
## Per-group equivalents of the whole-pop correlation_plot.R modes from the
## original monolithic pipeline. Each plot compares this group's own F_SNP
## (.het), F_ROH, D_STD/D_ROH (ROHRM vs GRM Inbreeding_Comparison), and
## ROH_sh (consensus intersect) — so the Trotter tab shows Trotter-only
## distributions, Pacer tab shows Pacer-only, and wholePop matches today.
RM_diag_rg="${canonical_dir}/Inbreeding_Comparison.${rg}.csv"
het_stats_rg="${het_rg_prefix}.het"
Froh_stats_rg="${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt"
conShare_rg="${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt"

## Strict GPA naming: wholePop outputs also carry the ".wholePop." marker.
froh_prefix="${OUTPUT_DIR}/divStats/coi_Froh_rmDiag_correlation.${rg}"
froh_cons_prefix="${OUTPUT_DIR}/divStats/coi_Froh_rmDiag_conShare_correlation.${rg}"
dbl_tag=".${rg}"

## --mode froh: COI vs F_ROH vs D_STD vs D_ROH
Rscript scripts/correlation_plot.R --mode froh "$RM_diag_rg" "$het_stats_rg" "$Froh_stats_rg" "$froh_prefix"
rclone -v copy "${froh_prefix}.pairplot.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## --mode froh-cons: adds ROH_sh via conShare
Rscript scripts/correlation_plot.R --mode froh-cons "$RM_diag_rg" "$het_stats_rg" "$Froh_stats_rg" "$conShare_rg" "$froh_cons_prefix"
rclone -v copy "${froh_cons_prefix}.pairplot.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## Focused doubleAnn plots: ROH_sh vs D_STD, F_ROH vs D_ROH, F_SNP vs D_ROH, F_SNP vs F_ROH.
## plot_correlation_withColorsAndShapes.R hard-codes its output filename to
## correlation_plot_<x>_vs_<y>_doubleAnn.png, so we rename after each call for
## Trotter / Pacer.
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$3;next}{print $0,a[$1]}' \
    <(cat "$conShare_rg" | sed 's/Percent_of_Consensus_ROH/ROH_sh/') \
    <(cat "$RM_diag_rg" | tr ',' '\t') \
    > "${OUTPUT_DIR}/divStats/rmdiag_conShare${dbl_tag}"
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' \
    "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize" \
    "${OUTPUT_DIR}/divStats/rmdiag_conShare${dbl_tag}" \
    > "${OUTPUT_DIR}/divStats/rmdiag_conShare${dbl_tag}_wBooksize"
Rscript "$scripts/plot_correlation_withColorsAndShapes.R" "${OUTPUT_DIR}/divStats/rmdiag_conShare${dbl_tag}_wBooksize" ROH_sh D_STD Phenotype Book_Size
    mv "${OUTPUT_DIR}/divStats/correlation_plot_ROH_sh_vs_D_STD_doubleAnn.png" "${OUTPUT_DIR}/divStats/correlation_plot_ROH_sh_vs_D_STD_doubleAnn${dbl_tag}.png"

rclone -v copy "${OUTPUT_DIR}/divStats/correlation_plot_ROH_sh_vs_D_STD_doubleAnn${dbl_tag}.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$5;next}{print $0,a[$1]}' "$Froh_stats_rg" <(cat "$RM_diag_rg" | tr ',' '\t') > "${OUTPUT_DIR}/divStats/rmdiag_Froh${dbl_tag}"
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize" "${OUTPUT_DIR}/divStats/rmdiag_Froh${dbl_tag}" > "${OUTPUT_DIR}/divStats/rmdiag_Froh${dbl_tag}_wBooksize"
Rscript "$scripts/plot_correlation_withColorsAndShapes.R" "${OUTPUT_DIR}/divStats/rmdiag_Froh${dbl_tag}_wBooksize" F_ROH D_ROH Phenotype Book_Size
    mv "${OUTPUT_DIR}/divStats/correlation_plot_F_ROH_vs_D_ROH_doubleAnn.png" "${OUTPUT_DIR}/divStats/correlation_plot_F_ROH_vs_D_ROH_doubleAnn${dbl_tag}.png"

rclone -v copy "${OUTPUT_DIR}/divStats/correlation_plot_F_ROH_vs_D_ROH_doubleAnn${dbl_tag}.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t"}NR==1{a[$2]="F_SNP";next}NR==FNR{a[$2]=$8;next}{print $0,a[$1]}' "$het_stats_rg" <(cat "$RM_diag_rg" | tr ',' '\t') > "${OUTPUT_DIR}/divStats/rmdiag_Fsnp${dbl_tag}"
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize" "${OUTPUT_DIR}/divStats/rmdiag_Fsnp${dbl_tag}" > "${OUTPUT_DIR}/divStats/rmdiag_Fsnp${dbl_tag}_wBooksize"
Rscript "$scripts/plot_correlation_withColorsAndShapes.R" "${OUTPUT_DIR}/divStats/rmdiag_Fsnp${dbl_tag}_wBooksize" F_SNP D_ROH Phenotype Book_Size
    mv "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_D_ROH_doubleAnn.png" "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_D_ROH_doubleAnn${dbl_tag}.png"

rclone -v copy "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_D_ROH_doubleAnn${dbl_tag}.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t"}NR==1{a[$2]="F_SNP";next}NR==FNR{a[$2]=$8;next}{print $0,a[$1]}' "$het_stats_rg" "${OUTPUT_DIR}/divStats/rmdiag_Froh${dbl_tag}" > "${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp${dbl_tag}"
awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$2]=$3;next}{if(a[$1])print $0,a[$1];}' "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize" "${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp${dbl_tag}" > "${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp${dbl_tag}_wBooksize"
Rscript "$scripts/plot_correlation_withColorsAndShapes.R" "${OUTPUT_DIR}/divStats/rmdiag_Froh_Fsnp${dbl_tag}_wBooksize" F_SNP F_ROH Phenotype Book_Size
    mv "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_F_ROH_doubleAnn.png" "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_F_ROH_doubleAnn${dbl_tag}.png"

rclone -v copy "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_F_ROH_doubleAnn${dbl_tag}.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
