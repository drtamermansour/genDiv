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

## The cross-group summaries currently live inline in genDiversity.sh's
## remaining Section 5/6 blocks (Froh_vs_ROHsh plots iterating over
## wholePop/twoGait/threeBooksize, merged_kin_sorted_top, etc.). They will
## move here in subsequent commits once the whole-pop KING / ROHRM / IBS /
## Euclidean work they depend on is itself moved out of the wrapper.

## For now this script is a no-op placeholder — it just sources common and
## logs its own entry so the wrapper can call it and the pipeline still runs
## end-to-end.
