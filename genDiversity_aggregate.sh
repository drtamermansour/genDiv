#!/usr/bin/env bash
# Cross-group aggregation stage for the genDiv pipeline.
# Runs once after genDiversity_per_group.sh has been invoked for all three
# groups. Contains only the blocks that truly require all three per-group
# runs: twoGait / threeBooksize per-sample concatenations, and the
# Froh-vs-ROHsh plots that consume them. Per-group wholePop-only analyses
# (F_ROH histograms, roh_high, gait/bookSize stratified summaries,
# merged_kin top-pair) live in genDiversity_per_group.sh under the
# if [[ "$rg" == "wholePop" ]] blocks.
set -eo pipefail
source "$(dirname "$0")/genDiversity_common.sh"

log "Cross-group aggregation stage"

roh_RG="${OUTPUT_DIR}/divStats/roh.L3"

##########################################
## Cross-group per-sample ROH_sh concatenations
##########################################
## twoGait   = Trotter + Pacer rows from each gait's own consensus-intersect file.
## threeBooksize = the six book-size subgroup rows (Trotter_LOW/MEDIUM/HIGH +
##                 Pacer_LOW/MEDIUM/HIGH).
## Each sample therefore appears exactly once in each concatenation, under the
## consensus of the gait (or book-size subgroup) it actually belongs to.
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
## table produced by per_group.sh wholePop (§8b); by the time aggregate.sh
## runs, that file already exists.
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
