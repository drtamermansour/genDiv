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
## Per-group sample lists (consumed by genDiversity_per_group.sh)
##########################################
## samples.${rg}.txt — 2-column FID\tIID files suitable for plink2 --keep or
## bcftools view -S (via `cut -f2` for a 1-column IID list when needed).
## wholePop = all samples in the post-QC LD-pruned fam (the authoritative set).
## Trotter / Pacer = intersection of the fam with gait labels. Samples without
## a gait label appear only in wholePop.
awk 'BEGIN{OFS="\t"}{print $1,$2}' "${pl1_pruned}.fam" > "${OUTPUT_DIR}/preprocess/samples.wholePop.txt"
for rg in Trotter Pacer; do
    awk -v rg="$rg" 'BEGIN{OFS="\t"} NR==FNR{fam[$2]=$1;next} fam[$2] && $3==rg {print fam[$2],$2}' \
        "${pl1_pruned}.fam" "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait" \
        > "${OUTPUT_DIR}/preprocess/samples.${rg}.txt"
done

## sample_groups.tsv — canonical IID→primary-group map for the reference set.
## One row per sample; wholePop membership is implicit (everyone is in wholePop).
## group ∈ {Trotter, Pacer, wholePop}. Samples with no gait label get group=wholePop.
{
    echo "# sample_groups.tsv — primary group per reference sample."
    echo "# One row per sample; wholePop membership is implicit."
    echo "# group ∈ {Trotter, Pacer, wholePop}. Samples with no gait label get group=wholePop."
    printf "IID\tgroup\n"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{gait[$2]=$3;next} {g=gait[$2]; if(g==""||g=="undefined") g="wholePop"; print $2,g}' \
        "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait" "${pl1_pruned}.fam"
} > "${OUTPUT_DIR}/preprocess/sample_groups.tsv"

log "Per-group sample lists: wholePop=$(wc -l < "${OUTPUT_DIR}/preprocess/samples.wholePop.txt") Trotter=$(wc -l < "${OUTPUT_DIR}/preprocess/samples.Trotter.txt") Pacer=$(wc -l < "${OUTPUT_DIR}/preprocess/samples.Pacer.txt")"



##########################################
## Effective autosomal genome length + autosomes.genome
##########################################
## calculate the effective autosomal genome length
awk '{print $1"\t"$4}' "$pl1_filtered".bim | grep "^chr" | grep -v "^chrX" > "$pl1_filtered".snp_pos.txt
aut_len=$(sort -k1,1 -k2,2n "$pl1_filtered".snp_pos.txt | \
        awk '{if ($1 == prev_chr) { gap = $2 - prev_pos; \
              if(gap > 0) {if (gap > 1000000) gap = 1000000; total += gap; }}\
              prev_chr=$1; prev_pos=$2} END {print total}') ## 2,261,547,402
echo $aut_len > ${OUTPUT_DIR}/divStats/effective_autosomal_genome_length.txt

## autosomes.genome — chromosome sizes table consumed by per_group.sh bedtools genomecov
awk '$1 ~ /^[0-9]+$/' "$reference_fai" | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2}' > "${OUTPUT_DIR}/divStats/autosomes.genome"

##########################################
## Whole-pop KING-robust kinship + IBS (consumed by per_group.sh Euclidean + correlation plots)
##########################################
group="gait"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --make-king-table 'counts' 'cols=+ibs1' \
    --output-chr 'chrM' --out "${OUTPUT_DIR}/divStats/filtered.LD_prune.king_$group"

kingkin="${OUTPUT_DIR}/divStats/filtered.LD_prune.king_$group.kin0"
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($10/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' <(tail -n+2 $kingkin) > ${kingkin%.kin0}.histo
rclone -v copy "$kingkin" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
rclone -v copy "${kingkin%.kin0}.histo" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## Likely first-degree relations (whole-pop filter; per-gait filters in per_group.sh)
head -n 1 "$kingkin" > "${OUTPUT_DIR}/divStats/related"
tail -n +2 "$kingkin" | sort -grk10,10 | awk '{if($10>0.177)print}' >> "${OUTPUT_DIR}/divStats/related"

## IBS augmentation of kin0: IBS1 = HET1_HOM2 + HET2_HOM1; IBS2 = N_SNPs - (HETHET + IBS0 + IBS1);
## IBS = (2*IBS2 + IBS1) / (2*N_SNPs).
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0,"IBS";next}{ibs1=$8+$9;ibs2=$5-($6+$7+ibs1);print $0,(2*ibs2+ibs1)/(2*$5)}' "$kingkin" > "${kingkin}.withIBS"

## Whole-pop KING kinship vs IBS correlation plot (one-shot)
Rscript "$scripts/plot_correlation.R" "${kingkin}.withIBS" KINSHIP IBS
rclone -v copy "${OUTPUT_DIR}/divStats/correlation_plot_KINSHIP_vs_IBS.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

##########################################
## A_e — Effective number of alleles (whole-pop diversity metric)
##########################################
## A_e = 1/Σ p_i^2 per SNP; reported as mean ± SD across SNPs, at three
## strata: the whole reference, per-gait (via --loop-cats 'PHENO1'), and per
## gait × book-size. This is a whole-population diversity summary, not a
## per-group reference consumed by GPA, so it lives here in shared.sh rather
## than per_group.sh.
########################################################
## 1. Effective number of alleles (\(A_{e}\)) 
## A_e represents the number of equally frequent alleles required to achieve the same level of expected heterozygosity (\(H_{e}\)) observed in a population
########################################################

## Formula
## \(A_{e} = \frac{1}{\sum p_{i}^{2}}\)
## where \(p_{i}\) is the frequency of the \(i^{th}\) allele

## Example Calculation
## Suppose a single locus has three alleles with the observed frequencies (0.6, 0.3, 0.1) in a population.
## 1. Calculate the squared frequencies: 0.36, 0.09, and 0.01
## 2. Calculate \(A_{e}\): 1/(0.36 + 0.09 + 0.01) = 1/0.46 = 2.17
## This result means that although there are 3 distinct alleles, the population's genetic diversity is equivalent to a population with only 2.17 equally frequent alleles. 

## Calculate \(A_{e}\) for each SNP
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --freq \
    --out "$pl1_pruned.freq_stats"
awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.afreq" > "$pl1_pruned.freq_stats.wholePop.afreq.Ae"
awk -v pop="wholePop" 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
    { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"pop, mean_Ae, "SD_Ae_in_"pop, sd_Ae } }' "$pl1_pruned.freq_stats.wholePop.afreq.Ae"
#Mean_Ae_in_wholePop     1.5616  SD_Ae_in_wholePop       0.321994

## Calculate \(A_{e}\) for each SNP per gait subpopulation
group="gait"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group \
    --loop-cats 'PHENO1' --freq \
    --out "$pl1_pruned.freq_stats"
#--loop-cats: Processing category 'Pacer' (271 samples).
#--loop-cats: Processing category 'Trotter' (271 samples).

## Calculate Mean \(A_{e}\) and standard deviation per gait subpopulation
for pop in Trotter Pacer; do
    awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.$pop.afreq" > "$pl1_pruned.freq_stats.$pop.afreq.Ae"
    awk -v pop=$pop 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
        { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"pop, mean_Ae, "SD_Ae_in_"pop, sd_Ae } }' "$pl1_pruned.freq_stats.$pop.afreq.Ae"
done
#Mean_Ae_in_Trotter      1.53063 SD_Ae_in_Trotter        0.339752
#Mean_Ae_in_Pacer        1.55017 SD_Ae_in_Pacer  0.331251

## Calculate \(A_{e}\) for each SNP per book size in each gait subpopulation
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize \
    --loop-cats 'PHENO1' --freq \
    --out "$pl1_pruned.freq_stats"
#--loop-cats: Processing category 'Pacer_HIGH' (75 samples).
#--loop-cats: Processing category 'Pacer_LOW' (90 samples).
#--loop-cats: Processing category 'Pacer_MEDIUM' (106 samples).
#--loop-cats: Processing category 'Trotter_HIGH' (58 samples).
#--loop-cats: Processing category 'Trotter_LOW' (96 samples).
#--loop-cats: Processing category 'Trotter_MEDIUM' (117 samples).

## Calculate Mean \(A_{e}\) and standard deviation per book size
for pop in Trotter Pacer; do
    for book in LOW MEDIUM HIGH;do
        awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.${pop}_${book}.afreq" > "$pl1_pruned.freq_stats.${pop}_${book}.afreq.Ae"
        awk -v gp=${pop}_${book} 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
            { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"gp, mean_Ae, "SD_Ae_in_"gp, sd_Ae } }' "$pl1_pruned.freq_stats.${pop}_${book}.afreq.Ae"
    done
done
#Mean_Ae_in_Pacer_LOW    1.55141 SD_Ae_in_Pacer_LOW      0.330248
#Mean_Ae_in_Pacer_MEDIUM 1.54804 SD_Ae_in_Pacer_MEDIUM   0.331937
#Mean_Ae_in_Pacer_HIGH   1.53805 SD_Ae_in_Pacer_HIGH     0.338842
#Mean_Ae_in_Trotter_LOW  1.53338 SD_Ae_in_Trotter_LOW    0.338901
#Mean_Ae_in_Trotter_MEDIUM       1.53176 SD_Ae_in_Trotter_MEDIUM 0.340108
#Mean_Ae_in_Trotter_HIGH 1.50661 SD_Ae_in_Trotter_HIGH   0.34827


Rscript scripts/effAllele_stats.R "${OUTPUT_DIR}" &> ${OUTPUT_DIR}/divStats/effAllele_stats.txt
Rscript scripts/plot_Ae.R "${OUTPUT_DIR}"
rclone -v copy ${OUTPUT_DIR}/divStats/effAllele_stats.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Ae/" --drive-shared-with-me
rclone -v copy ${OUTPUT_DIR}/divStats/Figure_Ae_BookSize.tiff "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Ae/" --drive-shared-with-me
#rclone -v copy ${OUTPUT_DIR}/divStats/Figure_Ae_BookSize.pdf "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Ae/" --drive-shared-with-me


##########################################
## FST between subpopulations (whole-pop differentiation metric)
##########################################
## `for group in sex gait bookSize; do plink2 --fst ...; done` computes
## autosomal FST for three groupings, then cats the per-grouping summaries
## into autosomal.fst.summary and feeds fst_stats.R for the adjusted-pairwise
## table. This is a one-shot whole-pop analysis — no per-group parameterisation.
##########################################
## 2. Fst between subpopulations (genders, gait types, and book sizes)
##########################################
## The fixation index can range from 0 to 1, where 0 means complete sharing of genetic material and 1 means no sharing. 
## For values equal to 1(meaning no sharing), scientists say that the populations are fixed.
## Effects of marker type and filtering criteria on QST-FST comparisons: https://pmc.ncbi.nlm.nih.gov/articles/PMC6894560/
for group in sex gait bookSize; do
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --pheno ${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.$group \
        --fst 'PHENO1' 'blocksize=2000' \
        --output-chr 'chrM' --out ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_$group
done

find ${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_*.summary -maxdepth 1 -type f | grep -v "\.x\." | xargs cat > ${OUTPUT_DIR}/divStats/autosomal.fst.summary
rclone -v copy ${OUTPUT_DIR}/divStats/autosomal.fst.summary "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me


Rscript scripts/fst_stats.R "${OUTPUT_DIR}" &> ${OUTPUT_DIR}/divStats/fst_stats.txt
rclone -v copy ${OUTPUT_DIR}/divStats/fst_stats.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me
