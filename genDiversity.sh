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
aut_len=$(cat "${OUTPUT_DIR}/divStats/effective_autosomal_genome_length.txt")

## Per-group stage. Each invocation produces reference files for one group.
## Currently a stub; subsequent commits move code here from the Section 5/6
## blocks below.
for rg in wholePop Trotter Pacer; do
    bash "$(dirname "$0")/genDiversity_per_group.sh" "$rg"
done

############################################
## x. Nucleotide diversity statistic (pi) -- This section is under development
############################################
## Nucleotide diversity is a population-level metric, the average number of differences between a pair of chromosomes, across all chromosome combinations within the population.
## This is distinct from simply measuring heterozygosity.
## Variant vs. Invariant Sites: Traditional pi calculations require knowledge of both variant and invariant sites (i.e., sequencing data).
##    With a SNP array, pi will be overestimated because it ignores the conserved (non-variable) parts of the genome (i.e., it is "SNP-based" diversity rather than a true "genomic" diversity.)



## Section 6 relatedness (whole-pop GRM + ROHRM cutoff loop + analysis_comparison
## + kinship-diff analysis) fully moved into per_group.sh §9. The wholePop
## invocation reproduces today's whole-pop outputs under .wholePop. suffixes.




## KING / IBS / PCA Euclidean / Euclidean-KING merge / cross-method correlation
## plot all moved to shared.sh (one-shot whole-pop pieces) and per_group.sh
## (per-group pieces). The merged_kin_sorted_top cross-reference moves to
## aggregate.sh next commit.

########################################################

## Cross-group aggregation (currently a no-op stub; future commits will
## move twoGait/threeBooksize summaries + Froh_vs_ROHsh plots + merged_kin_sorted_top
## here from the inline Section 5/6 blocks above).
bash "$(dirname "$0")/genDiversity_aggregate.sh"
