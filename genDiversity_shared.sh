#!/usr/bin/env bash
# Shared whole-pop preprocessing stage for the genDiv pipeline.
# Runs once. Contains Sections 1–3: data download, ID/sex cleanup, EquCab3
# remapping, SNP deduplication, QC filtering (MAF/HWE/missingness), and LD
# pruning. Outputs are consumed by genDiversity_per_group.sh and
# genDiversity_aggregate.sh.
set -eo pipefail
source "$(dirname "$0")/genDiversity_common.sh"

## Download genotyping data (PLINK: ped and map files)
log "Section 1: Downloading data"
module load rclone ## Loading rclone/1.65.1
mkdir -p "${OUTPUT_DIR}/SNPdata_iScan_Standardbred"
SNPdata="$(pwd)/${OUTPUT_DIR}/SNPdata_iScan_Standardbred"
#rclone lsd remote_UCDavis_GoogleDr: --drive-shared-with-me
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/SNP data - iScan_Standardbred" --drive-shared-with-me --include "USTA_Diversit*" $SNPdata/.

## Download metadata
mkdir -p "${OUTPUT_DIR}/Miscellaneous_documents_standardbred"
docs="$(pwd)/${OUTPUT_DIR}/Miscellaneous_documents_standardbred"
#rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/Miscellaneous documents_standardbred/USTA_Gait_BookSize_Assignments_Sex_Added.xlsx" --drive-shared-with-me $docs/.
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/updated_resources/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.xlsx" --drive-shared-with-me $docs/.
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/updated_resources/QC excluded samples/trotters_toExclude.lst" --drive-shared-with-me $docs/.
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/updated_resources/QC excluded samples/pacers_toExclude.lst" --drive-shared-with-me $docs/.

DOCS_DIR="$docs" python3 - <<'EOF'
import os
import pandas as pd
docs = os.environ["DOCS_DIR"]
df = pd.read_excel(os.path.join(docs, "USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.xlsx"), sheet_name="Sheet1")
df.to_csv(os.path.join(docs, "USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv"), index=False)
EOF

## create a tsv version of the metadata
cat $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv | tr ' ' '_' | tr ',' '\t' > $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.tsv

## QC and preprocessingx
mkdir -p ${OUTPUT_DIR}/preprocess
## check metadata
## Confrim that each sire show up in one book size
tail -n+2 $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv | \
    cut -d"," -f5,7 | sort -t"," -k2,2 | uniq > ${OUTPUT_DIR}/preprocess/sire_book
awk 'BEGIN{FS=","}{horses[$2]++}END{ \
	if(length(horses) < NR) { \
		print "Oops! We have these duplicate sires in the input BookSize file." > "/dev/stderr"; \
		for(h in horses) { if(horses[h]>1) print h } \
	} else { print "Good! No duplicate sires in the in the input BookSize file." > "/dev/stderr";} \
}' ${OUTPUT_DIR}/preprocess/sire_book | grep -Fwf - ${OUTPUT_DIR}/preprocess/sire_book || true

## Identify full siblings 
tail -n+2 $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv | \
    cut -d"," -f1-5,7-8 | sort -t"," -k3,3 -k6,6 -k7,7 | \
    awk 'BEGIN{FS=","}{ key = $6 OFS $7 } seen[key]++ { print prev_line ORS $0; next } { prev_line = $0 }' | cut -d, -f1 | paste - - > ${OUTPUT_DIR}/preprocess/full_siblings

## read genotypes (PLINK: bed + bim + fam files are writtin)
plink --file $SNPdata/USTA_Diversity_Study --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --output-chr 'M' --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex

## Update of ids (e.g., change HR15423_1 to HR15423)
cat ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.fam | tr ' ' '\t' | cut -f1-2 | grep "_" > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.cur_ids
sed 's/_.*//' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.cur_ids > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.new_ids
paste ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.cur_ids ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.new_ids > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.update_ids
plink --bfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --update-ids ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex.update_ids \
        --make-bed --output-chr 'M' --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex_updatedIDs

## Add sex metadata
cat ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex_updatedIDs.fam | tr ' ' '\t' | cut -f1-2 | tr '\t' ',' > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.ids
awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $1,$2,a[$2];}' \
    $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.ids > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.sex
plink --bfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study_noSex_updatedIDs --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --update-sex ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.sex \
        --make-bed --output-chr 'M' --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study

## Create ID lists
awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$1]=$3;next}{if(a[$2])print $1,$2,a[$2];}' \
    $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.ids | grep -vFwf <(cat $docs/{trotters,pacers}_toExclude.lst) > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ## 558

awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$1]=$5;next}{if(a[$2])print $1,$2,a[$2];}' \
    $docs/USTA_CuratedGait_BookSize_Assignments_with_Sires_and_Dams_CompositeBS.csv ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.ids > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize ## 576

awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$2]=$3;next}{if(a[$2])print $1,$2,a[$2]"_"$3}' \
    ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bookSize > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize ## 558

##########################################
## Remapping to EquCab3 coordinates
##########################################
log "Section 2: Remapping to EquCab3"
## The mapping file of Equine80select markers in EquCab3 coordinates
## chr \t pos \t snpID \t SNP_alleles \t genomic_alleles \t SNP_ref_alleles \t genomic_ref_allele \t allele_usage_decision
## This map allows updating the input alleles (using their IDs) into the SNP_alleles (i.e., Manifest alleles) or genomic_alleles (i.e., Positive Strand alleles = VCF alleles).
## In either case, their is a ref_allele to use in PLINK2
## check if the SNP alleles in BIM match those in the equCab3_map file
cat ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.bim | awk 'BEGIN{FS=OFS="\t"}{if($5 && $6){a[1]=$5;a[2]=$6;asort(a);print $2,a[1]","a[2]}}' > ${OUTPUT_DIR}/preprocess/tmpX_alleles_in_BIM.txt ## e.g., "UKUL1_ilmndup1  A,G"
tail -n+2 $equCab3_map | awk 'BEGIN{FS=OFS="\t"}{split($4, a, ",");asort(a);split($5, b, ",");asort(b);print $3,a[1]","a[2],b[1]","b[2]}' > ${OUTPUT_DIR}/preprocess/tmpX_alleles_in_MAP.txt ## e.g., "21962991_Curly_f_ilmndup1       A,G     A,G"
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$1])print $1,a[$1],$2,$3;}' \
    ${OUTPUT_DIR}/preprocess/tmpX_alleles_in_BIM.txt ${OUTPUT_DIR}/preprocess/tmpX_alleles_in_MAP.txt > ${OUTPUT_DIR}/preprocess/tmpX_compare_BIM_MAP.txt ## SNP_ID \t BIM_alleles \t MAP_SNP_alleles \t MAP_genomic_alleles
awk 'BEGIN{FS=OFS="\t"}{if($2!=$3)a+=1;if($2!=$4)b+=1;}END{print "mismatching SNP alleles:",a," mismatching genomic alleles:",b;}' ${OUTPUT_DIR}/preprocess/tmpX_compare_BIM_MAP.txt
## mismatching SNP alleles:        35961    mismatching genomic alleles:   36685
## Let us remove the ambiguous SNPs (A/T or C/G) from the analysis to avoid strand issues
awk 'BEGIN{FS=OFS="\t"}{if($4=="A,T" || $4=="T,A" || $4=="C,G" || $4=="G,C")print $3}' $equCab3_map > ${OUTPUT_DIR}/preprocess/ambiguous_snps.txt ## 262
## Also, let us remove the SNPs on unplaced Scaffolds
cat $equCab3_map | grep ^Un_NW | cut -f3 > ${OUTPUT_DIR}/preprocess/unplaced_snps.txt

## 1. select the variants to keep  
## 2. update chr/positions based on the equCab3_map
## 3. update -ve strand SNP alleles to postive strand version
tail -n+2 $equCab3_map | cut -f3 | grep -v -f <(cat ${OUTPUT_DIR}/preprocess/ambiguous_snps.txt ${OUTPUT_DIR}/preprocess/unplaced_snps.txt) > ${OUTPUT_DIR}/preprocess/snps_to_remap.txt ## 79314
tail -n+2 $equCab3_map | awk 'BEGIN{FS=OFS="\t"}{print $5}' | tr 'TCGA' 'AGCT' > ${OUTPUT_DIR}/preprocess/temp_pos_strand_complement.txt ## complementary genomic_alleles
paste <(tail -n+2 $equCab3_map) ${OUTPUT_DIR}/preprocess/temp_pos_strand_complement.txt | awk 'BEGIN{FS=OFS="\t"}{print $3,$9,$5}' | tr ',' '\t' > ${OUTPUT_DIR}/preprocess/pos_strand_alleles.txt ## SNP_ID \t complementary_genomic_alleles \t genomic_alleles
plink --bfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --extract ${OUTPUT_DIR}/preprocess/snps_to_remap.txt \
    --update-chr $equCab3_map 1 3 1 \
    --update-map $equCab3_map 2 3 1 \
    --update-alleles ${OUTPUT_DIR}/preprocess/pos_strand_alleles.txt \
    --make-bed --output-chr 'M' --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap ## input BIM has 79,259 ==> 76,841 remaining

## 4. update genomic alleles to fill in missing alleles (useless but just to be complete and make sure no snps will show up as mismtach in the next step)
tail -n+2 $equCab3_map | awk 'BEGIN{FS=OFS="\t"}{print $3,$5,$5}' | tr ',' '\t' > ${OUTPUT_DIR}/preprocess/genomic_alleles.txt ## SNP_ID \t genomic_alleles \t genomic_alleles
plink --bfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --update-alleles ${OUTPUT_DIR}/preprocess/genomic_alleles.txt \
    --make-bed --output-chr 'M' --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap

## check if the SNP alleles in BIM match genomic_alleles in the equCab3_map file after strand update
cat ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.bim | awk 'BEGIN{FS=OFS="\t"}{a[1]=$5;a[2]=$6;asort(a);print $2,a[1]","a[2]}' > ${OUTPUT_DIR}/preprocess/tmpX_alleles_in_remap.BIM.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$1])print $1,a[$1],$2,$3;}' \
    ${OUTPUT_DIR}/preprocess/tmpX_alleles_in_remap.BIM.txt ${OUTPUT_DIR}/preprocess/tmpX_alleles_in_MAP.txt > ${OUTPUT_DIR}/preprocess/tmpX_compare_remap.BIM_MAP.txt ## SNP_ID \t BIM_alleles \t MAP_SNP_alleles \t MAP_genomic_alleles
awk 'BEGIN{FS=OFS="\t";a=b=0;}{if($2!=$3)a+=1;if($2!=$4)b+=1;}END{print "mismatching SNP alleles:",a," mismatching genomic alleles:",b;}' ${OUTPUT_DIR}/preprocess/tmpX_compare_remap.BIM_MAP.txt
## mismatching SNP alleles:        35417    mismatching genomic alleles:   0

##########################################
## deduplication of SNPs based on chromosome and position
##########################################
mkdir -p ${OUTPUT_DIR}/dedup
# 1. compute per-SNP missingness:
plink --bfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --missing \
      --output-chr 'chrM' --out ${OUTPUT_DIR}/dedup/USTA_Diversity_Study.remap.missing
# 2. list positions that occur more than once (chr:bp repeated):
awk '{print $1":"$4}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.bim | sort | uniq -c | awk '$1>1{print $2}' > ${OUTPUT_DIR}/dedup/dup_positions.txt
# 3. extract SNP IDs at those duplicate positions:
# produce tab: SNP_ID <TAB> CHR:BP
awk 'BEGIN{FS=OFS="\t"} {print $2, $1":"$4}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.bim > ${OUTPUT_DIR}/dedup/bim_pos.tsv

# keep only rows where position is duplicated
grep -F -f ${OUTPUT_DIR}/dedup/dup_positions.txt ${OUTPUT_DIR}/dedup/bim_pos.tsv > ${OUTPUT_DIR}/dedup/duplicates_snps.tsv
# duplicates_snps.tsv: SNP_ID <TAB> CHR:BP (only positions that had >1 SNP)

# 4. join missingness to those SNPs and pick the best per position (lowest F_MISS = highest call rate)
# prepare a quick lookup of missingness: SNP_ID <TAB> F_MISS
awk 'BEGIN{OFS="\t"}NR>1{print $2, $5}' ${OUTPUT_DIR}/dedup/USTA_Diversity_Study.remap.missing.lmiss > ${OUTPUT_DIR}/dedup/snp_miss.tsv

# join: we want lines with SNP_ID, POS, F_MISS
awk 'BEGIN{FS=OFS="\t"} NR==FNR{miss[$1]=$2; next} {print $1,$2,miss[$1]}' ${OUTPUT_DIR}/dedup/snp_miss.tsv ${OUTPUT_DIR}/dedup/duplicates_snps.tsv > ${OUTPUT_DIR}/dedup/dup_with_miss.tsv
# dup_with_miss.tsv columns: SNP_ID  CHR:BP  F_MISS

# sort by position then by F_MISS ascending and pick the first SNP (best) per position
sort -k2,2 -k3,3n ${OUTPUT_DIR}/dedup/dup_with_miss.tsv | awk -F"\t" '{
  pos=$2;
  if(!(pos in seen)){ print $1"\t"$2"\t"$3; seen[pos]=1}
}' > ${OUTPUT_DIR}/dedup/best_per_pos.tsv
# best_per_pos.tsv: selected SNP_ID per duplicated position (the ones we keep)

# 5. produce a list of SNPs to remove (all duplicates except the selected ones):
# all duplicated SNP IDs:
cut -f1 ${OUTPUT_DIR}/dedup/duplicates_snps.tsv > ${OUTPUT_DIR}/dedup/all_dup_ids.txt
# selected to keep:
cut -f1 ${OUTPUT_DIR}/dedup/best_per_pos.tsv > ${OUTPUT_DIR}/dedup/keep_ids.txt
# produce remove list = setdiff(all_dup_ids - keep_ids)
grep -v -Fwf ${OUTPUT_DIR}/dedup/keep_ids.txt ${OUTPUT_DIR}/dedup/all_dup_ids.txt > ${OUTPUT_DIR}/preprocess/remove_dup_ids.txt

# 6. remove them with PLINK:
plink --bfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --exclude ${OUTPUT_DIR}/preprocess/remove_dup_ids.txt --make-bed \
    --output-chr 'chrM' --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.dedup
# 76841 variants loaded from .bim file.
# 71548 variants pass filters and QC.

## Convert PLINK.1 files to PLINK.2 binary format
plink2 --bfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.dedup --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --ref-allele 'force' $equCab3_map 7 3 1 --real-ref-alleles \
        --make-pgen --sort-vars \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2

##########################################
## Data exploration
##########################################
mkdir -p ${OUTPUT_DIR}/inspect
## --check-sex compares sex assignments in the input dataset with those imputed from chrX inbreeding coefficients 
## Preliminary run of --check-sex without removing PAR regions
plink2 --pfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
      --check-sex max-female-xf=$SEX_MAX_FEMALE_XF min-male-xf=$SEX_MIN_MALE_XF \
      --output-chr 'chrM' --out ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex
tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck | tr ' ' '\t' | cut -f3-5 | sort | uniq -c
#    285 1       1       OK
#    247 2       2       OK
#     44 2       NA      PROBLEM

## Histogram of X chromosome inbreeding coefficients (output of preliminary --check-sex)
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck) > ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck.histo

## Define Par regions and remove them (using high confidence males)
awk '$6 > 0.95 {print $1, $2}' ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck > ${OUTPUT_DIR}/inspect/hiConf_males.txt
#plink2 --pfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study --keep ${OUTPUT_DIR}/inspect/hiConf_males.txt --het --out ${OUTPUT_DIR}/inspect/male_het
cat ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.sex  | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"0"}' > ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.NoSex 
plink2 --pfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
       --keep ${OUTPUT_DIR}/inspect/hiConf_males.txt --update-sex ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.NoSex --geno-counts --out ${OUTPUT_DIR}/inspect/freqx_male
head -n1  ${OUTPUT_DIR}/inspect/freqx_male.gcount >  ${OUTPUT_DIR}/inspect/freqx_male.X.gcount
grep "^X" ${OUTPUT_DIR}/inspect/freqx_male.gcount >>  ${OUTPUT_DIR}/inspect/freqx_male.X.gcount
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($1=="chrX")a[$3]=$2;next}{$1=a[$2];print $0;}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.pvar ${OUTPUT_DIR}/inspect/freqx_male.X.gcount > ${OUTPUT_DIR}/inspect/freqx_male.X_ann.gcount
## Heterozygosity is seen until 2063653 which match our expectations (The PAB location on the X chromosome of EquCab2 is located at 1,175,430 bp)
## There is no detected coordinates for the tail PAR, thus I will use the position of the last marker + 1

## --check-sex with removal of PAR and noisy regions
awk -v par_end=$PAR_END_BP 'BEGIN{FS="\t"}{if($1<par_end || $6>5)print $2}' ${OUTPUT_DIR}/inspect/freqx_male.X_ann.gcount >  ${OUTPUT_DIR}/inspect/PAR_and_noise.list
plink2 --pfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
    --exclude ${OUTPUT_DIR}/inspect/PAR_and_noise.list \
    --check-sex max-female-xf=$SEX_MAX_FEMALE_XF min-male-xf=$SEX_MIN_MALE_XF \
    --output-chr 'chrM' --out ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex

tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck | tr ' ' '\t' | cut -f3-5 | sort | uniq -c
#    285 1       1       OK
#    247 2       2       OK
#     44 2       NA      PROBLEM
##  285 (male) and 247 (female) and 44 (NA). (Females are lost due to higher inbreeding)

## Histogram of X chromosome inbreeding coefficients (output of final --check-sex)
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck) > ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck.histo


## --het, --missing, --freq, --hardy
plink2 --pfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
      --het --missing --freq --hardy 'midp'  \
      --output-chr 'chrM' --out ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore

#--freq: Allele frequencies (founders only) written to ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.afreq .
#--missing: Sample missing data report written to ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.smiss .
#--missing: Variant missing data report written to ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.vmiss .
#--hardy midp: Autosomal Hardy-Weinberg report (founders only) written to ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy .
#--hardy midp: chrX Hardy-Weinberg report (founders only) written to ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy.x .
#Excluding 3581 variants on non-autosomes from --het.
#--het: done.
#Warning: 3937 variants skipped because they were monomorphic. 
#use --read-freq to provide more accurate allele frequency estimates.
#--het: Results written to ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het .

## heterozygosity (Postive estimates indicate high homozygosity while negative estimates indicate low homozygosity.)
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                END { for(i=bmin;i<=bmax;++i){if(i==0) print -1*size,size,a[i]/1;else if(i<0) print (i-1)*size,i*size,a[i]/1;else print i*size,(i+1)*size,a[i]/1 }}'  <(tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het) > ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het.histo

: <<'COMMENT'
-0.1    -0.08   3
-0.08   -0.06   6
-0.06   -0.04   17
-0.04   -0.02   26
-0.02   0.02    113
0.02    0.04    89
0.04    0.06    94
0.06    0.08    82
0.08    0.1     61
0.1     0.12    31
0.12    0.14    22
0.14    0.16    11
0.16    0.18    11
0.18    0.2     6
0.2     0.22    2
0.22    0.24    1
0.24    0.26    1
COMMENT

## Look as samples with heck-sex problem
paste ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het | cut -f1,2,4,6,9-12 | awk -F"\t" '{if($3=="NA")print}' | sort -t $'\t' -k4,4g > ${OUTPUT_DIR}/inspect/problem.sexcheck


## HWE
awk 'BEGIN{OFS="\t";}{ if($10<1e-50)a["1e-50 or less"]++;
                       else if($10<1e-40)a["1e-40:1e-50"]++; else if($10<1e-30)a["1e-30:1e-40"]++; \
                       else if($10<1e-20)a["1e-20:1e-30"]++; else if($10<1e-10)a["1e-10:1e-20"]++; \
                       else if($10<1e-9)a["1e-9:1e-10"]++; else if($10<1e-8)a["1e-8:1e-9"]++; \
                       else if($10<1e-7)a["1e-7:1e-8"]++; else if($10<1e-6)a["1e-6:1e-7"]++; \
                       else if($10<1e-5)a["1e-5:1e-6"]++; else if($10<1e-4)a["1e-4:1e-5"]++; \
                       else if($10<0.001)a["1e-3:1e-4"]++; else if($10<0.01)a["1e-2:1e-3"]++; \
                       else a["0.01 or more"]++; } \
                 END { for(i in a) print i,a[i] }'  <(tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy) | sort -g > ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy.histo

: <<'COMMENT'
1e-50 or less   17
1e-40:1e-50     11
1e-30:1e-40     25
1e-20:1e-30     79
1e-10:1e-20     592
1e-9:1e-10      197
1e-8:1e-9       259
1e-7:1e-8       419
1e-6:1e-7       575
1e-5:1e-6       713
1e-4:1e-5       1241
1e-3:1e-4       1987
1e-2:1e-3       4015
0.01 or more    57892
COMMENT

## Use this to further explore variants with extreme deviation from HWE:
cat ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy | awk '{if(NR==1)print}{if($10<1e-50)print}' >  ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy.lowHWE



## Variants with very low MAF
# Here are 2 different resolutions for a histogram of MAF:
awk -v size=0.01 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.afreq) > ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.frq.histo
awk -v size=0.001 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.afreq) > ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.frq.histo2

##########################################
## Final filtering based on missingness, MAF, and HWE
log "Section 3: Filtering (MAF=${MAF}, HWE=${HWE_PVAL}, missingness=${GENO_MISS})"
mkdir -p ${OUTPUT_DIR}/filtered

# Identify possible related dogs.
plink2 --pfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
       --king-cutoff $KING_CUTOFF \
       --out ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.1st_degree_relatives
# Make sure to include one of each known full siblings
tail -n+2 ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.1st_degree_relatives.king.cutoff.out.id | cut -f2 | \
    grep -vFwf - ${OUTPUT_DIR}/preprocess/full_siblings | cut -f1 | grep -Fwf - ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.psam > ${OUTPUT_DIR}/preprocess/full_siblings_rep
cat ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.1st_degree_relatives.king.cutoff.out.id ${OUTPUT_DIR}/preprocess/full_siblings_rep | cut -f1,2 > ${OUTPUT_DIR}/preprocess/relatives_toBeExcluded

# Run filtration
plink2 --pfile ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
      --remove ${OUTPUT_DIR}/preprocess/relatives_toBeExcluded \
      --hwe $HWE_PVAL 'midp' --geno $GENO_MISS --mind $GENO_MISS --maf $MAF --autosome \
      --real-ref-alleles --make-pgen --output-chr 'chrM' --out ${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered

#--remove: 560 samples remaining.
#560 samples (285 females, 275 males; 560 founders) remaining after main
#0 samples removed due to missing genotype data (--mind).
#--geno: 415 variants removed due to missing genotype data.
#--hwe midp: 2029 variants removed due to Hardy-Weinberg exact test (founders only).
# 7749 variants removed due to allele frequency threshold(s)
# 57829 variants remaining after main filters.

## Check final genotyping rate
plink2 --pfile ${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered --genotyping-rate \
       --out ${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered.genotyping_rate ## Total (hardcall) genotyping rate is 0.997906.

## Convert back to PLINK1 binary format for compatibility with other tools
pl1_filtered="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered"
plink2 --pfile ${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered \
      --real-ref-alleles --make-bed --output-chr 'chrM' --out $pl1_filtered

##########################################
## Convert to VCF format
plink2 --pfile ${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered \
      --real-ref-alleles --export vcf id-paste=iid --output-chr 'chrM' --out ${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered
vcf_filtered="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.vcf"


# check ref alleles and positions of the VCFs
grep -v '^##chrSet' $vcf_filtered | grep -E "^#|^chr" | bgzip --output $vcf_filtered.test.gz
tabix $vcf_filtered.test.gz
bcftools norm -c ws -f $ref $vcf_filtered.test.gz 1> $vcf_filtered.test.check.vcf 2> $vcf_filtered.test.check.log
## Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      57829/0/0/0/0/0/0
## REF/ALT total/modified/added:   57829/0/0

##########################################
## Convert to phased VCF format
## Prepare VCF for BEAGLE phasing and bcftools roh
grep -v '^##chrSet' $vcf_filtered | grep -E "^#|^chr" | grep -v "^chrX" | bgzip --output $vcf_filtered.auto.gz
tabix $vcf_filtered.auto.gz

## Double check the vcf file
bcftools norm \
  --rm-dup exact \
  -Oz \
  -o $vcf_filtered.norm.vcf.gz \
  $vcf_filtered.auto.gz ## Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      57829/0/0/0/0/0/0
tabix -p vcf $vcf_filtered.norm.vcf.gz

## Run BEAGLE (seed fixed for reproducibility)
beagle gt=$vcf_filtered.norm.vcf.gz out=$vcf_filtered.norm.phased nthreads=$nthreads seed=$BEAGLE_SEED
# Effective population size (Ne) is the number of individuals in an idealized population that would experience the same amount of genetic drift or inbreeding as the real, observed population. 
# we should provide this number as an input to Beagle when imputing few samples in the mating app.
grep "Estimated ne" $vcf_filtered.norm.phased.log | awk -F":" '{a+=$2}END{print "Ave. Estimated ne:",a/NR}' # Ave. Estimated ne: 2771.08
tabix -p vcf $vcf_filtered.norm.phased.vcf.gz

## Assess change in genotyping rate
plink2 --vcf $vcf_filtered.norm.vcf.gz --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --genotyping-rate --out $vcf_filtered.norm.genotyping_rate ## Total (hardcall) genotyping rate is 0.997906.
plink2 --vcf $vcf_filtered.norm.phased.vcf.gz --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --genotyping-rate --out $vcf_filtered.norm.phased.genotyping_rate ## Total (hardcall) genotyping rate is 1.

##########################################
## LD pruning to get independent variants for diversity calculations (& and output as PLINK1 binary format)
mkdir -p ${OUTPUT_DIR}/LD_pruned
plink2 --vcf $vcf_filtered.norm.phased.vcf.gz --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --indep-pairwise ${LD_WINDOW_KB}kb $LD_R2 \
       --real-ref-alleles --output-chr 'chrM' --out ${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.LD_lst ## 12147/57829 variants removed

pl1_pruned="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune"
plink2 --vcf $vcf_filtered.norm.phased.vcf.gz --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --psam ${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered.psam \
       --extract ${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.LD_lst.prune.in \
       --real-ref-alleles --make-bed --output-chr 'chrM' --out $pl1_pruned ## 45576 variants remaining


## Explore the LD-pruned dataset
plink2 --bfile $pl1_pruned --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --het --missing --freq --hardy 'midp'  \
      --output-chr 'chrM' --out ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune.explore

## check the change in (F) between: ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het ${OUTPUT_DIR}/inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune.explore.het | less ## F (i.e., measurement of inbreeding) decrease after pruning

## Check final genotyping rate of the LD-pruned dataset
plink2 --bfile $pl1_pruned --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --genotyping-rate --out $pl1_pruned.genotyping_rate ## Total (hardcall) genotyping rate is 0.997828.

## Convert to VCF format
plink2 --vcf $vcf_filtered.norm.phased.vcf.gz --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --extract ${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.LD_lst.prune.in \
       --real-ref-alleles --export vcf id-paste=iid --output-chr 'chrM' --out ${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.LD_prune 
vcf_pruned="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.norm.phased.LD_prune.vcf"

# check ref alleles and positions of the VCFs
grep -v '^##chrSet' $vcf_pruned | grep -E "^#|^chr" | bgzip --output $vcf_pruned.test.gz
tabix $vcf_pruned.test.gz
bcftools norm -c ws -f $ref $vcf_pruned.test.gz 1> $vcf_pruned.test.check.vcf 2> $vcf_pruned.test.check.log
## Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      45682/0/0/0/0/0/0
## REF/ALT total/modified/added:   45682/0/0


##########################################
## Whole-pop ROH calling, L1/L2/L3 filters, per-sample summaries
##########################################
##########################################
# 4E. ROH using bcftools/roh (Filtered dataset without LD pruning) -- This is the final approved approach
##########################################
## Run bcftools roh
bcftools roh -G30 --estimate-AF - $vcf_filtered.norm.phased.vcf.gz -o ${OUTPUT_DIR}/divStats/roh_out.txt
##Number of target samples: 560
##Number of --estimate-AF samples: 560
##Number of sites in the buffer/overlap: unlimited
##Number of lines total/processed: 57829/57829 (old: 58106/58106)
##Number of lines ${OUTPUT_DIR}/filtered/no AF/no alt/multiallelic/dup: 0/0/0/0/0

grep -E "^RG|^#" ${OUTPUT_DIR}/divStats/roh_out.txt > ${OUTPUT_DIR}/divStats/roh_out_RG.txt
## Summary stats by RG
awk 'BEGIN{print "IID\tNSEG\tKB\tKBAVG"} $1=="RG"{n[$2]++; sum[$2]+=$6} END{for (s in n) printf "%s\t%d\t%.2f\t%.2f\n", s, n[s], sum[s]/1000, (sum[s]/1000)/n[s]}' ${OUTPUT_DIR}/divStats/roh_out_RG.txt > ${OUTPUT_DIR}/divStats/roh_summary_by_RG.txt
awk 'NR > 1{ sum2 += $2; sum3 += $3; sum4 += $4 } END \
    { count = NR - 1; printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum2/count, sum3/count, sum4/count }' ${OUTPUT_DIR}/divStats/roh_summary_by_RG.txt
##Average Number of runs of homozygosity (NSEG) : 86.53
##Average of the total length of runs (kb) across all samples: 456,016.51
##Average of the average length of runs (KBAVG) across all samples: 5,256.67

## filtration to match the PLINK quality suggestions 
#Minimum ROH length (--homozyg-kb) 1000 kb
awk '/^#/ || $6 >= 1000000' ${OUTPUT_DIR}/divStats/roh_out_RG.txt > ${OUTPUT_DIR}/divStats/roh.L1.txt
#Minimum number of SNPs in ROH (--homozyg-snp) 50
awk '/^#/ || $7 >= 50' ${OUTPUT_DIR}/divStats/roh.L1.txt > ${OUTPUT_DIR}/divStats/roh.L2.txt
#Quality scores
awk -v size=2 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(grep -v "^#" ${OUTPUT_DIR}/divStats/roh.L2.txt) > ${OUTPUT_DIR}/divStats/roh.L2.histo 
awk '/^#/ || $8 >= 20' ${OUTPUT_DIR}/divStats/roh.L2.txt > ${OUTPUT_DIR}/divStats/roh.L3.txt

## Summary stats by RG after QC filtration && Stratify the file by the gait type
## The stats will be recalculated again later with Froh using "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt"
awk 'BEGIN{print "IID\tNSEG\tKB\tKBAVG"} $1=="RG"{n[$2]++; sum[$2]+=$6} END{for (s in n) printf "%s\t%d\t%.2f\t%.2f\n", s, n[s], sum[s]/1000, (sum[s]/1000)/n[s]}' ${OUTPUT_DIR}/divStats/roh.L3.txt > ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt
awk 'BEGIN{FS=OFS="\t";gait["IID"]="gait"}FNR==NR{gait[$2]=$3;next} {if(gait[$1])print $0,gait[$1];else print $0,"undefined";}' ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait ${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3.txt > ${OUTPUT_DIR}/divStats/roh.L3_gait.txt
INPUT_ROH="${OUTPUT_DIR}/divStats/roh.L3_gait.txt"
OUTPUT_FILE="${OUTPUT_DIR}/divStats/roh.L3_gait.sumStats.csv"
python scripts/summary_roh.py -i "$INPUT_ROH" -o "$OUTPUT_FILE"
rclone -v copy ${OUTPUT_DIR}/divStats/roh.L3_gait.sumStats.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
#Subgroup,          N,      NSEG,           KB,                         KBAVG
#Whole Population,  560,    54.58 +/- 9.43, 413086.61 +/- 103171.97,    7533.22 +/- 1210.63
#Pacer,             271,    50.34 +/- 7.11, 379860.14 +/- 81886.48,     7535.88 +/- 1224.94
#Trotter,           271,    59.28 +/- 9.02, 451519.43 +/- 105452.45,    7580.75 +/- 1166.36
#undefined,         18,     47.67 +/- 11.92,334702.03 +/- 138748.03,    6777.52 +/- 1454.36

##########################################
## Per-base consensus ROH (wholePop + gait + book-size subsets)
##########################################
## A per-base consensus ROH where ≥25% of samples are in ROH filtered by minimum size 500 kb and stratified by gait type
## With and without applying a smoothing function to the per-base coverage data to reduce noise before identifying consensus ROH regions
#roh_RG=${OUTPUT_DIR}/divStats/roh_out_RG
roh_RG=${OUTPUT_DIR}/divStats/roh.L3
# 1. Convert RG output → BED format
awk 'BEGIN{OFS="\t"} $1=="RG" {print $3, $4-1, $5, $2}' "${roh_RG}.txt" > "${roh_RG}.bed"
# 2. Ensure ROHs from the same sample do not double-count
cut -f4 "${roh_RG}.bed" | sort -u | while read S; do
  awk -v s="$S" '$4==s' "${roh_RG}.bed" | sort -k1,1 -k2,2n | bedtools merge -i - | awk -v s="$S" 'BEGIN{OFS="\t"}{print $1,$2,$3,s}'
done | sort -k1,1 -k2,2n > "${roh_RG}.merged_per_sample.wholePop.bed"

# subset the bed file for each subpopulation
for rg in "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do 
    grep "$rg" ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize | cut -f2 | grep -f - "${roh_RG}.merged_per_sample.wholePop.bed" > "${roh_RG}.merged_per_sample.${rg}.bed"
done

# 3. Calculate per-base ROH frequency (i.e., how many samples are in ROH at each base position)
awk '$1 ~ /^[0-9]+$/' $reference_fai | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2}' > ${OUTPUT_DIR}/divStats/autosomes.genome
for rg in "wholePop" "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do 
    bedtools genomecov -i "${roh_RG}.merged_per_sample.${rg}.bed" -g ${OUTPUT_DIR}/divStats/autosomes.genome -bg > "${roh_RG}.per_base_coverage.${rg}.bed"
    awk -v size=5 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($4/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                      END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  "${roh_RG}.per_base_coverage.${rg}.bed" > "${roh_RG}.per_base_coverage.${rg}.histo"
    # upload bed files
    # pause for now to save space
    #rclone -v copy ${roh_RG}.per_base_coverage.${rg}.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/freq/" --drive-shared-with-me

    # upload histo files
    rclone -v copy ${roh_RG}.per_base_coverage.${rg}.histo "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/freq/" --drive-shared-with-me
done

# 4. Identify consensus ROH regions (≥25% of samples in ROH) and merge adjacent regions (minimum size 500 kb)
# With and without appling a smoothing function which adjust the per-base coverage value of regions briding intervals with high coverage. The function would assign the average coverage of the region and two flanking regions to the bridged interval.
for rg in "wholePop" "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do
  num_samples=$(cut -f4 "${roh_RG}.merged_per_sample.${rg}.bed" | sort -u | wc -l)
  threshold=$(echo "$pct * $num_samples / 100" | bc -l)
  ## Find consensus before smoothing
  awk -v threshold=$threshold 'BEGIN{OFS="\t"} $4 >= threshold {print}' "${roh_RG}.per_base_coverage.${rg}.bed" > "${roh_RG}.consensus_${pct}pct.${rg}.bed"
  bedtools merge -i "${roh_RG}.consensus_${pct}pct.${rg}.bed" -c 4 -o mean | awk -v min_mb=$CONSENSUS_MIN_MB 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=min_mb)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${rg}.bed"

  # Summary stats of consensus ROH regions
  echo "==== consensus ROH in ≥${pct}% of ${rg} samples BEFORE SMOOTHING ======"
  awk -v rg="$rg" -v nsam="$num_samples" 'BEGIN{OFS=",";maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
    {print rg,"\nNo. of segments","Total length (KB)","Ave. length (KB)","Max % of samples in consensus","Average % of samples in consensus",\
    "\n"NR,sumLen,sumLen/NR,(maxConsen/nsam)*100"%",((sumSamples/NR)/nsam)*100"%"}' "${roh_RG}.consensus_${pct}pct.merged.${rg}.bed"

  # Smooth per-base coverage: only average a middle interval if it exactly bridges two adjacent intervals
  # and both flanking intervals are >= threshold while the middle < threshold.
  awk -v thr="$threshold" 'BEGIN{OFS="\t"} {chr[NR]=$1; st[NR]=$2; en[NR]=$3; cov[NR]=$4} END{
      for(i=1;i<=NR;i++){
          newcov=cov[i]
          if(i>1 && i<NR){
              # check perfect contiguity: prev_end == cur_start && cur_end == next_start
              if(en[i-1]==st[i] && en[i]==st[i+1]){
                  if(cov[i-1] >= thr && cov[i+1] >= thr && cov[i] < thr){
                      newcov = (cov[i-1] + cov[i] + cov[i+1]) / 3
                  }
              }
          }
          printf "%s\t%d\t%d\t%.6f\n", chr[i], st[i], en[i], newcov
      }
  }' "${roh_RG}.per_base_coverage.${rg}.bed" | awk 'BEGIN{OFS="\t"}{$4=$4+0;print}' > "${roh_RG}.per_base_coverage.${rg}.smoothed.bed"
  # Generate new consensus using the smoothed per-base coverage
  awk -v threshold="$threshold" 'BEGIN{OFS="\t"} $4 >= threshold {print}' "${roh_RG}.per_base_coverage.${rg}.smoothed.bed" > "${roh_RG}.consensus_${pct}pct.${rg}.smoothed.bed" ## recovered 72 more regions for wholePop
  bedtools merge -i "${roh_RG}.consensus_${pct}pct.${rg}.smoothed.bed" -c 4 -o mean | awk -v min_mb=$CONSENSUS_MIN_MB 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=min_mb)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed"
  
  # upload bed files
  # pause for now to save space
  #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${rg}.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
  #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

  # Summary stats of consensus ROH regions
  echo "==== consensus smoothed ROH in ≥${pct}% of ${rg} samples AFTER SMOOTHING ======"
  awk -v rg="$rg" -v nsam="$num_samples" 'BEGIN{OFS=",";maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
   {print rg,"\nNo. of segments","Total length (KB)","Ave. length (KB)","Max % of samples in consensus","Average % of samples in consensus",\
    "\n"NR,sumLen,sumLen/NR,(maxConsen/nsam)*100"%",((sumSamples/NR)/nsam)*100"%"}' "${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed"
 
done > ${roh_RG}.consensus_${pct}pct.summary.txt
rclone -v copy ${roh_RG}.consensus_${pct}pct.summary.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
grep -A3 "AFTER SMOOTHING" ${roh_RG}.consensus_${pct}pct.summary.txt

## intersect ROH regions of each sample aganist the consensus ROH regions (in each "rg")
## Output: the census length and percentage in each sample (file for each "rg")
roh_RG="${OUTPUT_DIR}/divStats/roh.L3"
for rg in "wholePop" "Trotter" "Pacer" "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH"; do
    consensus_bed=${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed
    consensus_size=$(awk 'BEGIN{sum=0} {sum+=($3-$2)} END {print sum}' ${consensus_bed})
    bed_perSample="${roh_RG}.merged_per_sample.${rg}.bed"   ## no need to loop on ${rg} here. It should be the same if you always used "wholePop" 
    echo -e "IID\tTotal_ROH_in_Consensus_region(bp)\tPercent_of_Consensus_ROH" > ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt
    cut -f4 "${bed_perSample}" | sort -u | while read S; do
      awk -v s="$S" '$4==s' "${bed_perSample}" | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b "${consensus_bed}" | awk -v s="$S" -v cs="$consensus_size" 'BEGIN{OFS="\t"}{size+=($3-$2)} END {print s, size, (size/cs)*100}'
    done >> ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt
    #rclone -v copy ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
done

## merge the intersection (ROH_share) with Trotter/Pacer consensus
## ouput: the census length and percentage in each sample aganist its own gait consensus
head -n1  ${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt > ${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt
for rg in "Trotter" "Pacer";do 
    tail -n+2 ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt;
done >> ${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt  

## merge the intersection (ROH_share) with Trotter_booksize/Pacer_booksize consensus
## ouput: the census length and percentage in each sample aganist its own gait_bookSize consensus
head -n1  ${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt > ${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt
for rg in "Trotter_LOW" "Trotter_MEDIUM" "Trotter_HIGH" "Pacer_LOW" "Pacer_MEDIUM" "Pacer_HIGH";do 
    tail -n+2 ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt;
done >> ${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt  

##########################################
## Effective autosomal genome length
##########################################
## calculate the effective autosomal genome length
awk '{print $1"\t"$4}' "$pl1_filtered".bim | grep "^chr" | grep -v "^chrX" > "$pl1_filtered".snp_pos.txt
aut_len=$(sort -k1,1 -k2,2n "$pl1_filtered".snp_pos.txt | \
        awk '{if ($1 == prev_chr) { gap = $2 - prev_pos; \
              if(gap > 0) {if (gap > 1000000) gap = 1000000; total += gap; }}\
              prev_chr=$1; prev_pos=$2} END {print total}') ## 2,261,547,402
echo $aut_len > ${OUTPUT_DIR}/divStats/effective_autosomal_genome_length.txt
