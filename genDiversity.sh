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
     ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het > ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait

INPUT_HET="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait.sumStats.csv"
python scripts/summary_het.py -i "$INPUT_HET" -o "$OUTPUT_FILE"

## generate a summary table of heterozygosity and inbreeding coefficient in the three book size in the two subpopulations and the whole cohort
awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$2]=$3;next}{if(a[$2])print $0,a[$2];else print $0,"undefined";}' \
     ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het > ${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize

INPUT_HET="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het.wGait_bookSize.sumStats.csv"
python scripts/summary_het.py -i "$INPUT_HET" -o "$OUTPUT_FILE"

rclone -v copy ${OUTPUT_DIR}/divStats --include "filtered.LD_prune.het_stats.het*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/het_and_COI/" --drive-shared-with-me



##########################################
# 4. Runs of homozygosity (ROH)
##########################################
group="gait"
case="Trotter"

##########################################
## 4A. ROH using Plink (Pruned dataset) -- This is not actual useful. Only for testing the effect of pruning and using fewer count of markers 
##########################################
## Plink --homozyg
## By default, only runs of homozygosity containing at least 100 SNPs, and of total length ≥ 1000 kilobases, are noted. You can change these minimums with --homozyg-snp and --homozyg-kb, respectively.
## By default, a ROH must have at least one SNP per 50 kb on average; change this bound with --homozyg-density.
## By default, if two consecutive SNPs are more than 1000 kb apart, they cannot be in the same ROH; change this bound with --homozyg-gap.
## By default, a ROH can contain an unlimited number of heterozygous calls; you can impose a limit with --homozyg-het. (This flag was silently ignored by PLINK 1.07.)
## By default, the scanning window contains 50 SNPs; change this with --homozyg-window-snp.
## By default, a scanning window hit can contain at most 1 heterozygous call and 5 missing calls; change these limits with --homozyg-window-het and --homozyg-window-missing, respectively.
## By default, for a SNP to be eligible for inclusion in a ROH, the hit rate of all scanning windows containing the SNP must be at least 0.05; change this threshold with --homozyg-window-threshold.
## Due to how the scanning algorithm works, it is possible for a reported run of homozygosity to be adjacent to a few unincluded homozygous variants. This is generally harmless, but if you wish to extend the ROH to include them, use the 'extend' modifier. (Note that the --homozyg-density bound can prevent extension, and --homozyg-gap affects which variants are considered adjacent.)

plink --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group $case \
        --homozyg 'extend' \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group


## Outputs of --homozyg include: 
## .hom (run-of-homozygosity list)
##  FID, IID, PHE (henotype value), CHR, SNP1, SNP2, POS1, POS2, KB, NSNP, DENSITY, PHOM (% of homozygous), PHET (% of heterozygous), PMISS (% of missing), SENSITIVITY, SPECIFICITY

## .hom.indiv (sample-based runs-of-homozygosity report) which has the following columns: 
## FID	Family ID
## IID	Within-family ID
## PHE	Phenotype value
## NSEG	Number of runs of homozygosity
## KB	Total length of runs (kb)
## KBAVG	Average length of runs (kb)
## Calculate summary stats from .hom.indiv 
awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/count, sum5/count, sum6/count }' ${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 15.63 
##Average of the total length of runs (kb) across all samples: 164,427.23
##Average of the average length of runs (KBAVG) across all samples: 10,439.76


## Rscript that plots the correlation between "KB" and "KBAVG" from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.roh_$group.hom"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_pruned/" --drive-shared-with-me

##########################################
## 4B. ROH using Plink (Filtered dataset without LD pruning)
##########################################
plink --bfile "$pl1_filtered" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group $case \
        --homozyg 'extend' \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group

awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/count, sum5/count, sum6/count }' ${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 33.35
##Average of the total length of runs (kb) across all samples: 345,017.70
##Average of the average length of runs (KBAVG) across all samples: 10,305.44


## Rscript that plots the correlation between  KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_$group.hom"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_${OUTPUT_DIR}/filtered/" --drive-shared-with-me

##########################################
## 4C. ROH using Plink (Filtered dataset without LD pruning (with group option))
##########################################
plink --bfile "$pl1_filtered" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group $case \
        --homozyg group 'extend' \
        --homozyg-window-snp 20 \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group

awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/count, sum5/count, sum6/count }' ${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 35.46
##Average of the total length of runs (kb) across all samples: 364,115.99
##Average of the average length of runs (KBAVG) across all samples: 10,211.80

## Rscript that plots the correlation between  KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.group_roh_$group.hom"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_filtered_gp/" --drive-shared-with-me

############################################
## 4D. Studying ROH using Howard et al. (2016) approach -- This was for testing only and is not currently in use
############################################
mkdir -p $work_dir/${OUTPUT_DIR}/rep_ROHRM
rohrm_dir="$work_dir/${OUTPUT_DIR}/rep_ROHRM"

## ROH Analysis with Sub-populations and Phenotypes
## Having a headerless tab-separated file with 3 columns: The 2nd column has subject ids matching the VCF and the 3rd column has the a binary phenotype, let us do the following:
## 1. ROH Calling (Per Individual): We will implement a scanner that checks hap1 == hap2. Any contiguous stretch of matching haplotypes longer than the cutoff (e.g., 1 Mb) is flagged as an ROH.
## 2. Island Detection: We will calculate the frequency of ROHs at every SNP, find the "Top 5%" cutoff, and merge contiguous high-frequency SNPs into "Islands".
##    i.e., an "ROH Island" is defined strictly as a contiguous block of SNPs where every single SNP is in the Top 5% of frequencies.
## 3. Phenotype Integration: the analysis will be repeated for sub-populations defined in the phenotype file).
## 4. make a plot to show the ROH frequency across the genome for the whole population and each sub-population.
## 5. Calc the average (±SD) proportion of the genome in a ROH for the whole population and each sub-population.

roh_mb_cutoff=1.0  # in Megabases (Mb)
phenotypes="${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait"
python $scripts/ROH_analysis.py $vcf_filtered.norm.phased.vcf.gz $phenotypes $roh_mb_cutoff "$rohrm_dir" > $rohrm_dir/ROH_analysis.$roh_mb_cutoff.log

rclone -v copy $rohrm_dir/ROH_Frequency_Plot.$roh_mb_cutoff.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
rclone -v copy $rohrm_dir/ROH_Islands_Detailed.$roh_mb_cutoff.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
#rclone -v copy $rohrm_dir/ROH_Subpop_Stats.$roh_mb_cutoff.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
rclone -v copy $rohrm_dir/ROH_analysis.$roh_mb_cutoff.log "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me


tail -n+2 $rohrm_dir/ROH_Islands_Detailed.1.0.csv | awk -F"," '{sum += $5; array[NR] = $5} END {
    mean = sum / NR;
    for (x=1; x<=NR; x++) {
        sumsq += ((array[x] - mean)**2);
    }
    std_dev = sqrt(sumsq / NR); # Population standard deviation
    # For sample standard deviation, use sqrt(sumsq / (NR - 1)) if NR > 1
    min_snp_cutoff = mean - (2 * std_dev)

    print "Mean: " mean;
    print "Standard Deviation: " std_dev
    print "2 SD below Mean: " min_snp_cutoff
}' ## Mean: 28.6454 // Standard Deviation: 41.1817 // 2 SD below Mean: -53.7179 //There is no need to filter based on this criterion as it results in a negative value. As an alternative, we can use 3 as a minimum SNP count threshold for defining ROH islands.

## Filter ROH islands with at least 3 SNPs
awk -F"," 'NR==1 || $5 >= 3' $rohrm_dir/ROH_Islands_Detailed.1.0.csv > $rohrm_dir/Filtered_ROH_Islands_Detailed.1.0.csv
## for each sub-population in column 1, calculate the number of ROH regions (NR), the total length of ROH islands (sum of $4-$3), the average length, and the maximum and average SNP frequency (column 6)
awk -F"," 'NR>1 && $1!="" {
    grp=$1
    if (!(grp in seen)) { seen[grp]=1; order[++norder]=grp }
    len = ($4 - $3) + 0
    freq = ($6 + 0)
    count[grp]++
    sumlen[grp] += len
    sumfreq[grp] += freq
    if (!(grp in maxfreq) || freq > maxfreq[grp]) maxfreq[grp] = freq
}
END {
    print "Group,NR,Total_len_bp,Avg_len_bp,Max_freq,Avg_freq"
    for(i=1;i<=norder;i++) {
        g = order[i]
        printf "%s,%d,%.0f,%.2f,%.6f,%.6f\n", g, count[g], sumlen[g], sumlen[g]/count[g], maxfreq[g], sumfreq[g]/count[g]
    }
}' $rohrm_dir/Filtered_ROH_Islands_Detailed.1.0.csv > $rohrm_dir/Filtered_ROH_Subpop_Stats.csv

rclone -v copy $rohrm_dir/Filtered_ROH_Subpop_Stats.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/Howard_reimp/" --drive-shared-with-me
#python $scripts/ROH_analysis_withFilter.py filtered.norm.phased.vcf.gz phenotypes.txt $roh_mb_cutoff > Filtered_ROH_analysis.$roh_mb_cutoff.log



## Rscript that plots the correlation between KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
## Similar analysis will be done later after calculation of related matrices
roh_indiv="${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt" ## to read KB and KBAVG
het_stats="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_summary_by_RG_L3"
Rscript scripts/correlation_plot.R --mode basic $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me


############################################
## 5. F_ROH statistic (currently calculated based on bacftools roh)
############################################
## F_ROH is an inbreeding coefficient based on runs of homozygosity
## Standard practice is to calculate F_{ROH} statistics on all valid ROHs (>1Mb), while restricting "Islands" (signatures of selection) to only the most robust regions.
## per-sample F_ROH = (sum length of ROH for that individual) / (total autosomal genome length).


awk -v aut_len=$aut_len 'BEGIN{FS=OFS="\t";}NR==1{print $0,"F_ROH";next} {print $0, ($3*1000)/aut_len}' ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt > ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.txt

## Froh-vs-ROHsh plots moved to genDiversity_aggregate.sh; they consume the
## twoGait / threeBooksize concatenations produced there.

############################################
## x. Nucleotide diversity statistic (pi) -- This section is under development
############################################
## Nucleotide diversity is a population-level metric, the average number of differences between a pair of chromosomes, across all chromosome combinations within the population.
## This is distinct from simply measuring heterozygosity. 
## Variant vs. Invariant Sites: Traditional pi calculations require knowledge of both variant and invariant sites (i.e., sequencing data). 
##    With a SNP array, pi will be overestimated because it ignores the conserved (non-variable) parts of the genome (i.e., it is "SNP-based" diversity rather than a true "genomic" diversity.)



############################################
## 6. Relatedness work
############################################
log "Section 6: Relatedness analysis"
mkdir -p $work_dir/${OUTPUT_DIR}/rep_ROHRM
rohrm_dir="$work_dir/${OUTPUT_DIR}/rep_ROHRM"
group="gait"
phenotypes="${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait"


############################################
## Genomic Relatedness using Plink
############################################
## Create a standard GRM matrix 
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --make-rel square \
    --output-chr 'chrM' --out $rohrm_dir/filtered.LD_prune.GRM_$group

############################################
## ROH-based Relatedness (REPLICATE C++ code of Howard et al.)
############################################
## The C++ code you provided reveals a very specific logic: it relies on phased genotypes (distinguishing between Maternal vs. Paternal alleles), 
## and it calculates relationships based on exact haplotype matching within dynamically defined genomic windows.
# Stage 1: Loads the VCF, converts positions to Megabases (Mb), and splits genotypes into hap1/hap2.
#   To properly implement this, we previously converted the plink files into VCF, select autosomes, bcftools norm to remove duplicate markers, and phase by Beagle.
# Stage 2 (The Geometry Layer): Replicates the logic of the "ROH_Index class" in C++ to identify all valid genomic windows based on physical distance (Mb). 
#   It iterates through every SNP i and finds the window of SNPs that fits within the ROH Cutoff (e.g., 1 Mb). (Sliding Window with a step of 1 SNP)
    # Note: We are not checking for ROH. We are defining the windows only to be used in the next stages.
    # Note2: this implementation is faster than C++ which uses nested loops. Here, we use a single loop with a moving "end" pointer. This reduces the complexity from O(N^2) to O(N).
# Stage 3 (Window Filtering; The Statistical Layer): Discard windows that are too "sparse" (i.e., has too few SNPs) based on the statistical distribution.
    # a. Calculate the Mean and Standard Deviation (SD) of the "Number of SNPs" across all windows.
    # b. Define a threshold: Cutoff = Mean - (User_Threshold * SD).
    # c. Discard any window where NumSNP < Cutoff.
# Stage 4 (The Kernel): Calculate the relationship matrix. This is the hardest part to optimize in Python.
    # a. Extract the vector of genotypes for that window.
    # b. Expand genotypes into haplotypes (Paternal/Maternal).
    # c. Check similarity: If either haplotypes from animal i == either haplotypes from animal j, score 1.0 for any match; otherwise 0.0.
    #.   Note: when $i = j$ (The Diagonal): This function scores the inbreeding coefficient based on ROH: if hap1 == hap2, score 1.0; else 0.0. 
    # d. Accumulate these scores into the global matrix.
    # Note: Again, we optimize the C++ nested loops by using NumPy Broadcasting instead.
# Stage 5 (Normalization): Scale the matrix and save it.
: <<'COMMENT'
Feature,            Standard GRM (VanRaden),                    ROH GRM (Howard et al.)
Input Data,         "Unphased Genotypes (0, 1, 2)",             Phased Haplotypes (0|0, 0|1, 1|0, 1|1)
Resolution,         Single Nucleotide (SNP),                    Multi-Megabase Window (Window)
Matching Logic,     Allele Sharing (IBS),                       Exact String Matching (IBD)
Sensitivity,        Very tolerant of mutation/recombination.,   Very strict. One mismatch breaks the link.
Biological Signal,  Captures Deep/Ancient Relatedness.,         Captures Recent Relatedness.
Diagonal,           Heterozygosity-based Inbreeding.,           ROH-based Inbreeding (FROH​).
COMMENT

############################################
## Compare Genomic Relatedness and ROH-based Relatedness for each window size
## The script outputs "Robust_Matrix_Comparison_Enhanced.png", "Pairwise_Differences.csv", and "Inbreeding_Comparison.csv"
############################################
roh_threshold=$ROH_THRESHOLD_SD
for roh_mb_cutoff in $ROH_CUTOFFS; do
    roh_sd_label="${roh_threshold%.*}SD"
    subfolder="roh_${roh_mb_cutoff%.*}Mb.Threshold_${roh_sd_label}"
    log "Running ROHRM: ${roh_mb_cutoff} Mb cutoff, threshold ${roh_sd_label}"

    run_python $scripts/ROHRM_Creator.py $vcf_filtered.norm.phased.vcf.gz $roh_mb_cutoff $roh_threshold "$rohrm_dir" \
        > $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold.log

    run_python $scripts/analysis_comparison.py \
        $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold \
        $rohrm_dir/filtered.LD_prune.GRM_$group $phenotypes "$rohrm_dir"

    # center column 7 (Difference) around the column's mean and save to new column centered_Kinship_diff
    awk -F, 'BEGIN{OFS=FS=","} NR==1{hdr=$0; next} {sum+=$7; n++; lines[n]=$0; vals[n]=$7} \
             END{ mean = (n?sum/n:0); print hdr, "centered_Kinship_diff"; \
                  for(i=1;i<=n;i++) printf "%s%s%.8f\n", lines[i], OFS, vals[i]-mean }' \
        $rohrm_dir/Pairwise_Differences.csv > $rohrm_dir/Pairwise_Differences.csv.tmp \
        && mv $rohrm_dir/Pairwise_Differences.csv.tmp $rohrm_dir/Pairwise_Differences.csv

    # generate histograms of kinship distributions
    awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                        END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_Std.histo
    awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                        END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_ROH.histo
    awk -v size=0.01 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                        END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' \
        <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_diff.histo

    # extra pair analysis for the primary ROH cutoff only
    if [[ "$roh_mb_cutoff" == "$PRIMARY_ROH_MB" ]]; then
        ## animal pairs with high positive kinship difference (ROH-based kinship > standard kinship)
        awk 'BEGIN{FS=","} NR==1{print;next}{if($8>0.1) print}' \
            $rohrm_dir/Pairwise_Differences.csv > ${OUTPUT_DIR}/divStats/high_positive_kinship_diff.csv
        ## animal pairs with high negative kinship difference (standard kinship > ROH-based kinship)
        awk 'BEGIN{FS=","} NR==1{print;next}{if($8<-0.15) print}' \
            $rohrm_dir/Pairwise_Differences.csv > ${OUTPUT_DIR}/divStats/high_negative_kinship_diff.csv
        ## average centered kinship difference between and within gait groups
        awk 'BEGIN{FS=","} /Trotter/ && /Pacer/ {sum+=$8; n++} END{print "Ave diff Trotter-Pacer:", sum/n}' $rohrm_dir/Pairwise_Differences.csv
        awk 'BEGIN{FS=","} /Trotter/ && !/Pacer/ {sum+=$8; n++} END{print "Ave diff Trotter-Trotter:", sum/n}' $rohrm_dir/Pairwise_Differences.csv
        awk 'BEGIN{FS=","} !/Trotter/ && /Pacer/ {sum+=$8; n++} END{print "Ave diff Pacer-Pacer:", sum/n}' $rohrm_dir/Pairwise_Differences.csv
    fi

    # upload results
    upload $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png "Relatedness/$subfolder/"
    upload $rohrm_dir/Pairwise_Differences.csv "Relatedness/$subfolder/"
    upload $rohrm_dir/Inbreeding_Comparison.csv "Relatedness/$subfolder/"

    # move results to a per-cutoff subfolder
    mkdir -p $rohrm_dir/$subfolder
    mv $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png $rohrm_dir/Inbreeding_Comparison.csv \
       $rohrm_dir/Pairwise_Differences.* $rohrm_dir/$subfolder/
done




## KING / IBS / PCA Euclidean / Euclidean-KING merge / cross-method correlation
## plot all moved to shared.sh (one-shot whole-pop pieces) and per_group.sh
## (per-group pieces). The merged_kin_sorted_top cross-reference moves to
## aggregate.sh next commit.

########################################################

## Cross-group aggregation (currently a no-op stub; future commits will
## move twoGait/threeBooksize summaries + Froh_vs_ROHsh plots + merged_kin_sorted_top
## here from the inline Section 5/6 blocks above).
bash "$(dirname "$0")/genDiversity_aggregate.sh"
