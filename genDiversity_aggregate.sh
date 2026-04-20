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

## Additional cross-group outputs still live inline in genDiversity.sh:
## Froh-vs-ROHsh plots (for gp in wholePop twoGait threeBooksize) and the
## merged_kin_sorted_top across ROHRM/KING/PCA. They'll move here in a later
## commit once the whole-pop KING / IBS / Euclidean code they depend on is
## itself moved out of the wrapper.
