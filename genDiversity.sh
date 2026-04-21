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

############## Stats on diversity ##################
log "Section 5: Diversity statistics"
mkdir -p ${OUTPUT_DIR}/divStats



##########################################
## 3. Expected and observed heterozygosity and inbreeding coefficient
##########################################
## An inbreeding coefficient (COI) is a measure of the probability that an individual will have two copies of an allele that are identical by descent from a common ancestor.
## A higher COI means more predictability of traits but also a greater risk of genetic health problems due to inbreeding depression

## generate a summary table of heterozygosity and inbreeding coefficient in the two subpopulations and the whole cohort
awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' \
     ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.wholePop.het > ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait

INPUT_HET="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait.sumStats.csv"
python scripts/summary_het.py -i "$INPUT_HET" -o "$OUTPUT_FILE"

## generate a summary table of heterozygosity and inbreeding coefficient in the three book size in the two subpopulations and the whole cohort
awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' \
     ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.wholePop.het > ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize

INPUT_HET="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize.sumStats.csv"
python scripts/summary_het.py -i "$INPUT_HET" -o "$OUTPUT_FILE"

rclone -v copy ${OUTPUT_DIR}/divStats --include "filtered.LD_prune.het_stats.wholePop.het*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/het_and_COI/" --drive-shared-with-me


## Rscript that plots the correlation between KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
## Similar analysis will be done later after calculation of related matrices
roh_indiv="${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.wholePop.txt" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.wholePop.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_summary_by_RG_L3"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me


## F_ROH statistic: whole-pop F_ROH summary is produced by per_group.sh wholePop
## (roh_summary_by_RG_L3_Froh.wholePop.txt). Histograms, roh_high.csv, and
## gait/book-size stratified summaries + Froh-vs-ROHsh plots live in
## genDiversity_aggregate.sh.


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
