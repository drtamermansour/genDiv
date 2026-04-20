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
samples_rg="${OUTPUT_DIR}/preprocess/samples.${rg}.txt"
aut_len=$(cat "${OUTPUT_DIR}/divStats/effective_autosomal_genome_length.txt")

mkdir -p "${OUTPUT_DIR}/divStats"

## Subsequent commits move per-group code here:
##   5b — PCA (currently 3 hardcoded blocks in genDiversity.sh:~29–110)
##   5c — PCA + book-size overlay, FST book-size-within-gait, KING filter,
##        Euclidean distance, cross-method correlation plots
##   5d — NEW GPA-proposal metrics: per-group afreq, het, F_ROH, GRM, ROHRM,
##        analysis_comparison
