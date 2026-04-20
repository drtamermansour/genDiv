#!/usr/bin/env bash
# Cross-group aggregation stage for the genDiv pipeline.
# Runs once after genDiversity_per_group.sh has been invoked for all three
# groups. Consumes per-group outputs and produces summaries / plots that
# compare or combine them (twoGait / threeBooksize views, F_ROH-vs-ROHsh
# plots, merged kinship summary, etc.).
set -eo pipefail
source "$(dirname "$0")/genDiversity_common.sh"

log "Cross-group aggregation stage"

## Canonical paths produced by genDiversity_shared.sh and genDiversity_per_group.sh.
pl1_pruned="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune"
vcf_filtered="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.vcf"
phenotypes="${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait"
roh_RG="${OUTPUT_DIR}/divStats/roh.L3"

##########################################
## Cross-group per-sample ROH_sh concatenations
##########################################
## twoGait   = Trotter + Pacer rows from each gait's own consensus-intersect file.
## threeBooksize = the six book-size subgroup rows (Trotter_LOW/MEDIUM/HIGH +
##                 Pacer_LOW/MEDIUM/HIGH).
## Each sample therefore appears exactly once in each concatenation, under the
## consensus of the gait (or book-size subgroup) it actually belongs to. These
## files feed the Froh-vs-ROHsh plots in the inline Section-5 block of
## genDiversity.sh (iterating `for gp in wholePop twoGait threeBooksize`).
head -n1 "${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt" \
    > "${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt"
for rg in Trotter Pacer; do
    tail -n+2 "${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt"
done >> "${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt"

head -n1 "${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt" \
    > "${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt"
for rg in Trotter_LOW Trotter_MEDIUM Trotter_HIGH Pacer_LOW Pacer_MEDIUM Pacer_HIGH; do
    tail -n+2 "${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt"
done >> "${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt"

##########################################
## F_ROH vs ROH_sh plots (consume twoGait/threeBooksize summaries built above)
##########################################
## One plot per aggregation view. $froh is the gait_bookSize-stratified F_ROH
## table produced inline in genDiversity.sh's Section-5 block; by the time
## aggregate.sh runs, that file already exists.
froh="${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt"
for gp in wholePop twoGait threeBooksize; do
    conShare="${roh_RG}.perSample_intersect_${gp}_consensus_${pct}pct.summary.txt"
    output_file="${OUTPUT_DIR}/divStats/Froh_vs_ROHsh_${gp}.png"
    python scripts/roh_plot.py "$conShare" "$froh" "$output_file"
    rclone -v copy "$output_file" --drive-shared-with-me "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/"
    output_prefix="${OUTPUT_DIR}/divStats/normalized_ROHsh_${gp}"
    run_python scripts/roh_histograms.py --metric ratio "$conShare" "$froh" "$output_prefix"
    upload "$output_prefix".histogram.png "Froh/"
    upload "$output_prefix".density.png   "Froh/"
    output_prefix2="${OUTPUT_DIR}/divStats/ROHshared_${gp}"
    run_python scripts/roh_histograms.py --metric shared "$conShare" "$froh" "$output_prefix2"
    upload "$output_prefix2".histogram.png "Froh/"
done &> "${OUTPUT_DIR}/divStats/roh_sh.log"

## Still inline in genDiversity.sh: merged_kin_sorted_top across ROHRM/KING/PCA.
## That one depends on the whole-pop KING / Pairwise_Differences produced in
## the still-inline Section-6 block; it'll move here when that block is
## extracted.
