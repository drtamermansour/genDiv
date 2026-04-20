#!/usr/bin/env bash
# Per-group metrics stage for the genDiv pipeline.
# Invoked once per $rg ∈ {wholePop, Trotter, Pacer} by genDiversity.sh. Produces
# per-group reference files consumed by GPA's downstream report pipeline.
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

## Group-specific naming helpers. Keep wholePop's bare filenames unchanged
## during Step 5 so existing downstream consumers don't break; Step 7 will
## rename them to .wholePop. per the strict GPA convention.
if [[ "$rg" == "wholePop" ]]; then
    rg_tag=""          # inserted before ".eigenvec" etc. — empty for wholePop
else
    rg_tag=".${rg}"
fi

mkdir -p "${OUTPUT_DIR}/divStats"

##########################################
## PCA
##########################################
pca_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune${rg_tag}.pca"
if [[ "$rg" == "wholePop" ]]; then
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
           --real-ref-alleles --autosome --pca 'allele-wts' \
           --output-chr 'chrM' --out "$pca_prefix"
    n_pcs=6
else
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
           --keep "$samples_rg" --autosome --pca \
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
if [[ "$rg" == "wholePop" ]]; then
    out_png="${OUTPUT_DIR}/divStats/pca_plot_BookSize.png"
else
    out_png="${OUTPUT_DIR}/divStats/pca_plot_BookSize.${rg}.png"
fi
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

## Subsequent commits move more per-group code here:
##   5c — PCA + COI overlay (depends on per-group het file from 5d),
##        FST book-size-within-gait, KING filter, Euclidean distance,
##        cross-method correlation plots
##   5d — NEW GPA-proposal metrics: per-group afreq, het, F_ROH, GRM, ROHRM,
##        analysis_comparison

##########################################
## FST of book-size subgroups within this gait (Trotter / Pacer only)
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
## PCA: Color samples by COI (whole-pop het file from shared.sh)
##########################################
# Step 5d will switch each rg to its own --read-freq-derived .het; for now all
# three groups overlay COI from the whole-pop filtered.LD_prune.het_stats.het.
het_stats_whole="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"
awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 "$het_stats_whole") "$pca_prefix.eigenvec" > "$pca_prefix.eigenvec.wCOI"
eigenvec_suffix="wCOI"; color_column="COI"
if [[ "$rg" == "wholePop" ]]; then
    out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.png"
else
    out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.${rg}.png"
fi
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" "$n_pcs" numeric
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
