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

## Whole-pop KING outputs produced by shared.sh; per_group.sh consumes these
## for the per-gait related filter (§3.5 below) and for the Euclidean + KING
## merge + cross-method correlation plot (§10–§11).
group="gait"
kingkin="${OUTPUT_DIR}/divStats/filtered.LD_prune.king_${group}.kin0"
kingkin_wIBS="${kingkin}.withIBS"

## Group-specific naming helper. Keep wholePop's bare filenames unchanged so
## downstream consumers don't break; Step 7 will unify on the strict GPA
## convention (always ".wholePop." inserted) once GPA is ready.
if [[ "$rg" == "wholePop" ]]; then
    rg_tag=""
else
    rg_tag=".${rg}"
fi

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
## 4. Heterozygosity (F_SNP) using group-specific AF via --read-freq
##########################################
het_rg_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.${rg}"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$samples_rg" \
    --het 'cols=fid,hom,het,nobs,f' \
    --read-freq "$pruned_afreq" \
    --output-chr 'chrM' --out "$het_rg_prefix"
rclone -v copy "${het_rg_prefix}.het" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/het_and_COI/" --drive-shared-with-me

##########################################
## 5. PCA COI overlay (uses this group's own .het)
##########################################
awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 "${het_rg_prefix}.het") "$pca_prefix.eigenvec" > "$pca_prefix.eigenvec.wCOI"
eigenvec_suffix="wCOI"; color_column="COI"
if [[ "$rg" == "wholePop" ]]; then
    out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.png"
else
    out_png="${OUTPUT_DIR}/divStats/pca_plot_inbreeding.${rg}.png"
fi
Rscript scripts/pca_plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png" "$n_pcs" numeric
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

##########################################
## 6. bcftools roh + L1/L2/L3 filter (group-subset VCF so AF is group-specific)
##########################################
## For wholePop we use the whole-pop phased VCF directly; for Trotter / Pacer
## we create a group-subset phased VCF that is ALSO reused by ROHRM in §9.
if [[ "$rg" == "wholePop" ]]; then
    group_vcf="${vcf_filtered}.norm.phased.vcf.gz"
else
    group_vcf="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.${rg}.vcf.gz"
    bcftools view -S <(cut -f2 "$samples_rg") --force-samples \
        "${vcf_filtered}.norm.phased.vcf.gz" -Oz -o "$group_vcf"
    bcftools index -t "$group_vcf"
fi

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
## 7. Per-base consensus ROH (this group + its book-size subgroups)
##########################################
## Build a per-sample merged bed from the group's own L3 ROH calls, then for
## each "sub" of this group compute per-base ROH coverage, apply the ≥pct
## threshold, merge/smooth, and intersect per-sample with the consensus.
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
    # Book-size subsets reuse this gait's per-sample bed, filtered to book-size membership
    for sub in "${rg}_LOW" "${rg}_MEDIUM" "${rg}_HIGH"; do
        grep "$sub" "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize" | cut -f2 \
            | grep -f - "${roh_RG}.merged_per_sample.${rg}.bed" \
            > "${roh_RG}.merged_per_sample.${sub}.bed"
    done
fi

# 7.3 Per-base coverage + histograms + consensus threshold + smoothing + per-sample intersect
# autosomes.genome produced once in shared.sh; reuse.
{
for sub in "${subs[@]}"; do
    bed_perSample="${roh_RG}.merged_per_sample.${sub}.bed"
    bedtools genomecov -i "$bed_perSample" -g "${OUTPUT_DIR}/divStats/autosomes.genome" -bg > "${roh_RG}.per_base_coverage.${sub}.bed"
    awk -v size=5 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($4/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                          END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' "${roh_RG}.per_base_coverage.${sub}.bed" > "${roh_RG}.per_base_coverage.${sub}.histo"
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
done

##########################################
## 8. F_ROH summary (appends F_ROH column using shared effective genome length)
##########################################
awk -v aut_len="$aut_len" 'BEGIN{FS=OFS="\t";}NR==1{print $0,"F_ROH";next} {print $0, ($3*1000)/aut_len}' \
    "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.${rg}.txt" > "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt"
rclone -v copy "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/" --drive-shared-with-me

##########################################
## 9. GRM + ROHRM + analysis_comparison → D_SNP, G_SNP, D_ROH, G_ROH
##########################################
## Per-group standard GRM (vanRaden) via plink2 --make-rel + per-group ROHRM
## (Howard re-implementation). analysis_comparison.py writes
## Inbreeding_Comparison.csv (IID, D_STD, D_ROH, Phenotype) and
## Pairwise_Differences.csv (ID1, ID2, Pheno1, Pheno2, Kinship_Std, Kinship_ROH,
## Difference). Column schemas match the current whole-pop outputs so GPA's
## fixed-column-index reads keep working.
grp_workdir="${OUTPUT_DIR}/rep_ROHRM/perGroup_${rg}"
mkdir -p "$grp_workdir" "${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD"

grm_prefix="${grp_workdir}/filtered.LD_prune.GRM.${rg}"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$samples_rg" \
    --make-rel 'square' \
    --read-freq "$pruned_afreq" \
    --output-chr 'chrM' --out "$grm_prefix"

rohrm_mb="$PRIMARY_ROH_MB"
rohrm_sd="$ROH_THRESHOLD_SD"
rohrm_prefix="${grp_workdir}/ROHRM.rohMinSize_${rohrm_mb}.rohThreshold_${rohrm_sd}"
run_python "$scripts/ROHRM_Creator.py" "$group_vcf" "$rohrm_mb" "$rohrm_sd" "$grp_workdir" \
    > "${rohrm_prefix}.log"

run_python "$scripts/analysis_comparison.py" \
    "$rohrm_prefix" \
    "$grm_prefix" \
    "$phenotypes" "$grp_workdir"

canonical_dir="${OUTPUT_DIR}/rep_ROHRM/roh_1Mb.Threshold_3SD"
mv "${grp_workdir}/Inbreeding_Comparison.csv" "${canonical_dir}/Inbreeding_Comparison.${rg}.csv"
mv "${grp_workdir}/Pairwise_Differences.csv"  "${canonical_dir}/Pairwise_Differences.${rg}.csv"
rclone -v copy "${canonical_dir}/Inbreeding_Comparison.${rg}.csv" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
rclone -v copy "${canonical_dir}/Pairwise_Differences.${rg}.csv"  "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me

##########################################
## 10. Per-gait "related" filter (Trotter / Pacer only)
##########################################
## shared.sh's "related" file is the whole-pop list of first-degree pairs
## (KING kinship > 0.177). Here we keep just the rows where a gait-labeled
## sample appears. For wholePop the shared.sh file already covers the entire
## population, so no subsetting needed.
if [[ "$rg" != "wholePop" ]]; then
    grep "$rg" "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.${group}" | cut -f2 \
        | grep -Fwf - "${OUTPUT_DIR}/divStats/related" \
        > "${OUTPUT_DIR}/divStats/related_${rg}"
fi
