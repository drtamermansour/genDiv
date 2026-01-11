## update the conda package manager itself
conda update conda

## update all packages in an environment
conda update --all

## clean up unused packages and caches
conda clean --all

## create a new conda environment for genomic related GWAS tools
mamba create -n grGWAS
conda activate grGWAS
mamba install conda-forge::openpyxl conda-forge::pandas
mamba install -c bioconda plink plink2 bcftools gcta bedtools
mamba install -c conda-forge r-base=4.5.2 r-ggplot2=4.0.1 r-gridextra=2.3 r-qqman=0.1.9 r-viridis=0.6.5 r-reshape2=1.4.5 r-ggally=2.4.0

## Create the working directory of the project
mkdir -p $HOME/genDiv && cd $HOME/genDiv
mkdir -p scripts
scripts="$HOME/genDiv/scripts"

## Download genotyping data
module load rclone ## Loading rclone/1.65.1
mkdir -p SNPdata_iScan_Standardbred
SNPdata="$(pwd)/SNPdata_iScan_Standardbred"
rclone lsd remote_UCDavis_GoogleDr: --drive-shared-with-me
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/SNP data - iScan_Standardbred" --drive-shared-with-me --include "USTA_Diversit*" $SNPdata/.

## Download metaadata
mkdir -p Miscellaneous_documents_standardbred
docs="$(pwd)/Miscellaneous_documents_standardbred"
rclone -v copy "remote_UCDavis_GoogleDr:STR_Imputation_2025/Miscellaneous documents_standardbred/USTA_Gait_BookSize_Assignments_Sex_Added.xlsx" --drive-shared-with-me $docs/.

python3 - <<'EOF'
import pandas as pd
df = pd.read_excel("Miscellaneous_documents_standardbred/USTA_Gait_BookSize_Assignments_Sex_Added.xlsx", sheet_name="Sheet1")
df.to_csv("Miscellaneous_documents_standardbred/USTA_Gait_BookSize_Assignments_Sex_Added.csv", index=False)
EOF

## QC and preprocessing
mkdir -p preprocess
plink --file $SNPdata/USTA_Diversity_Study --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --output-chr 'M' --out preprocess/USTA_Diversity_Study_noSex

## Update of ids (e.g., change HR15423_1 to HR15423)
cat preprocess/USTA_Diversity_Study_noSex.fam | tr ' ' '\t' | cut -f1-2 | grep "_" > preprocess/USTA_Diversity_Study_noSex.cur_ids
sed 's/_.*//' preprocess/USTA_Diversity_Study_noSex.cur_ids > preprocess/USTA_Diversity_Study_noSex.new_ids
paste preprocess/USTA_Diversity_Study_noSex.cur_ids preprocess/USTA_Diversity_Study_noSex.new_ids > preprocess/USTA_Diversity_Study_noSex.update_ids
plink --bfile preprocess/USTA_Diversity_Study_noSex --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --update-ids preprocess/USTA_Diversity_Study_noSex.update_ids \
        --make-bed --output-chr 'M' --out preprocess/USTA_Diversity_Study_noSex_updatedIDs

## Add sex metadata
cat preprocess/USTA_Diversity_Study_noSex_updatedIDs.fam | tr ' ' '\t' | cut -f1-2 | tr '\t' ',' > preprocess/USTA_Diversity_Study.ids
awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $1,$2,a[$2];}' \
    $docs/USTA_Gait_BookSize_Assignments_Sex_Added.csv preprocess/USTA_Diversity_Study.ids > preprocess/USTA_Diversity_Study.sex
plink --bfile preprocess/USTA_Diversity_Study_noSex_updatedIDs --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --update-sex preprocess/USTA_Diversity_Study.sex \
        --make-bed --output-chr 'M' --out preprocess/USTA_Diversity_Study

## Create ID lists
awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$1]=$3;next}{if(a[$2])print $1,$2,a[$2];}' \
    $docs/USTA_Gait_BookSize_Assignments_Sex_Added.csv preprocess/USTA_Diversity_Study.ids > preprocess/USTA_Diversity_Study.gait

awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$1]=$5;next}{if(a[$2])print $1,$2,a[$2];}' \
    $docs/USTA_Gait_BookSize_Assignments_Sex_Added.csv preprocess/USTA_Diversity_Study.ids > preprocess/USTA_Diversity_Study.bookSize

##########################################
## New remapping to EquCab3 coordinates
##########################################
## The mapping file of Equine80select markers in EquCab3 coordinates 
## chr \t pos \t snpID \t SNP_alleles \t genomic_alleles \t SNP_ref_alleles \t genomic_ref_allele \t allele_usage_decision
## This map allows remapping to the same SNP alleles or PosStrand_alleles (i.e., VCF alleles). 
## In either case, their is a ref_allele to use in PLINK2
equCab3_map=$(pwd)/../Equine80select_remapper/matchingSNPs_binary_consistantMapping.EquCab3_map
## check if the SNP alleles in BIM match those in the equCab3_map file
cat preprocess/USTA_Diversity_Study.bim | awk 'BEGIN{FS=OFS="\t"}{if($5 && $6){a[1]=$5;a[2]=$6;asort(a);print $2,a[1]","a[2]}}' > preprocess/tmpX_alleles_in_BIM.txt
cat $equCab3_map | awk 'BEGIN{FS=OFS="\t"}{split($4, a, ",");asort(a);split($5, b, ",");asort(b);print $3,a[1]","a[2],b[1]","b[2]}' > preprocess/tmpX_alleles_in_MAP.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$1])print $1,a[$1],$2,$3;}' \
    preprocess/tmpX_alleles_in_BIM.txt preprocess/tmpX_alleles_in_MAP.txt > preprocess/tmpX_compare_BIM_MAP.txt ## SNP_ID \t BIM_alleles \t MAP_SNP_alleles \t MAP_genomic_alleles
awk 'BEGIN{FS=OFS="\t"}{if($2!=$3)a+=1;if($2!=$4)b+=1;}END{print "mismatching SNP alleles:",a," mismatching genomic alleles:",b;}' preprocess/tmpX_compare_BIM_MAP.txt
## mismatching SNP alleles:        36180    mismatching genomic alleles:   36900
## Let us remove the ambiguous SNPs (A/T or C/G) from the analysis to avoid strand issues
awk 'BEGIN{FS=OFS="\t"}{if($4=="A,T" || $4=="T,A" || $4=="C,G" || $4=="G,C")print $3}' $equCab3_map > preprocess/ambiguous_snps.txt ## 279

## 1. select the variants to keep  
## 2. update chr/positions based on the equCab3_map
## 3. update -ve strand SNP alleles to postive strand version
cut -f3 $equCab3_map | grep -v -f preprocess/ambiguous_snps.txt > preprocess/snps_to_remap.txt
awk 'BEGIN{FS=OFS="\t"}{print $5}' $equCab3_map | tr 'TCGA' 'AGCT' > preprocess/temp_pos_strand_complement.txt ## complementary genomic_alleles
paste $equCab3_map preprocess/temp_pos_strand_complement.txt | awk 'BEGIN{FS=OFS="\t"}{print $3,$9,$5}' | tr ',' '\t' > preprocess/pos_strand_alleles.txt ## SNP_ID \t complementary_genomic_alleles \t genomic_alleles
plink --bfile preprocess/USTA_Diversity_Study --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --extract preprocess/snps_to_remap.txt \
    --update-chr $equCab3_map 1 3 \
    --update-map $equCab3_map 2 3 \
    --update-alleles preprocess/pos_strand_alleles.txt \
    --make-bed --output-chr 'M' --out preprocess/USTA_Diversity_Study.remap ## 79259 ==> 77359 remaining (39458 variants update)

## 4. update genomic alleles to fill in missing alleles (useless but just to be complete and make sure no snps will show up as mismtach in the next step)
cat $equCab3_map | awk 'BEGIN{FS=OFS="\t"}{print $3,$5,$5}' | tr ',' '\t' > preprocess/genomic_alleles.txt ## SNP_ID \t genomic_alleles \t genomic_alleles
plink --bfile preprocess/USTA_Diversity_Study.remap --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --update-alleles preprocess/genomic_alleles.txt \
    --make-bed --output-chr 'M' --out preprocess/USTA_Diversity_Study.remap ## 79259 ==> 77359 remaining (39458 variants update)

## check if the SNP alleles in BIM match genomic_alleles in the equCab3_map file after strand update
cat preprocess/USTA_Diversity_Study.remap.bim | awk 'BEGIN{FS=OFS="\t"}{a[1]=$5;a[2]=$6;asort(a);print $2,a[1]","a[2]}' > preprocess/tmpX_alleles_in_remap.BIM.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$1])print $1,a[$1],$2,$3;}' \
    preprocess/tmpX_alleles_in_remap.BIM.txt preprocess/tmpX_alleles_in_MAP.txt > preprocess/tmpX_compare_remap.BIM_MAP.txt ## SNP_ID \t BIM_alleles \t MAP_SNP_alleles \t MAP_genomic_alleles
awk 'BEGIN{FS=OFS="\t";a=b=0;}{if($2!=$3)a+=1;if($2!=$4)b+=1;}END{print "mismatching SNP alleles:",a," mismatching genomic alleles:",b;}' preprocess/tmpX_compare_remap.BIM_MAP.txt
## mismatching SNP alleles:        35683    mismatching genomic alleles:   0

##XXXXXXXXX 4. set ref alleles
#awk 'BEGIN{FS=OFS="\t"}{if($4!=$5)print $3,$4,$5}' $equCab3_map | tr ',' '\t' > preprocess/SNP_alleles.txt ## snpID \t SNP_alleles \t genomic_alleles
plink --bfile preprocess/USTA_Diversity_Study.remap --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --update-alleles preprocess/SNP_alleles.txt \
    --a2-allele $equCab3_map 7 3 --real-ref-alleles \
    --make-bed \
    --output-chr 'M' --out preprocess/USTA_Diversity_Study.remap.refAlleles ## --a2-allele: 77630 assignments made.

##########################################
## deduplication of SNPs based on chromosome and position
##########################################
mkdir -p dedup
# 1. compute per-SNP missingness:
plink --bfile preprocess/USTA_Diversity_Study.remap --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --missing \
      --output-chr 'chrM' --out dedup/USTA_Diversity_Study.remap.missing
# 2. list positions that occur more than once (chr:bp repeated):
awk '{print $1":"$4}' preprocess/USTA_Diversity_Study.remap.bim | sort | uniq -c | awk '$1>1{print $2}' > dedup/dup_positions.txt
# 3. extract SNP IDs at those duplicate positions:
# produce tab: SNP_ID <TAB> CHR:BP
awk 'BEGIN{FS=OFS="\t"} {print $2, $1":"$4}' preprocess/USTA_Diversity_Study.remap.bim > dedup/bim_pos.tsv

# keep only rows where position is duplicated
grep -F -f dedup/dup_positions.txt dedup/bim_pos.tsv > dedup/duplicates_snps.tsv
# duplicates_snps.tsv: SNP_ID <TAB> CHR:BP (only positions that had >1 SNP)

# 4. join missingness to those SNPs and pick the best per position (lowest F_MISS = highest call rate)
# prepare a quick lookup of missingness: SNP_ID <TAB> F_MISS
awk 'BEGIN{OFS="\t"}NR>1{print $2, $5}' dedup/USTA_Diversity_Study.remap.missing.lmiss > dedup/snp_miss.tsv

# join: we want lines with SNP_ID, POS, F_MISS
awk 'BEGIN{FS=OFS="\t"} NR==FNR{miss[$1]=$2; next} {print $1,$2,miss[$1]}' dedup/snp_miss.tsv dedup/duplicates_snps.tsv > dedup/dup_with_miss.tsv
# dup_with_miss.tsv columns: SNP_ID  CHR:BP  F_MISS

# sort by position then by F_MISS ascending and pick the first SNP (best) per position
sort -k2,2 -k3,3n dedup/dup_with_miss.tsv | awk -F"\t" '{
  pos=$2;
  if(!(pos in seen)){ print $1"\t"$2"\t"$3; seen[pos]=1}
}' > dedup/best_per_pos.tsv
# best_per_pos.tsv: selected SNP_ID per duplicated position (the ones we keep)

# 5. produce a list of SNPs to remove (all duplicates except the selected ones):
# all duplicated SNP IDs:
cut -f1 dedup/duplicates_snps.tsv > dedup/all_dup_ids.txt
# selected to keep:
cut -f1 dedup/best_per_pos.tsv > dedup/keep_ids.txt
# produce remove list = setdiff(all_dup_ids - keep_ids)
grep -v -Fwf dedup/keep_ids.txt dedup/all_dup_ids.txt > preprocess/remove_dup_ids.txt

# 6. remove them with PLINK:
plink --bfile preprocess/USTA_Diversity_Study.remap --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --exclude preprocess/remove_dup_ids.txt --make-bed \
    --output-chr 'chrM' --out preprocess/USTA_Diversity_Study.remap.dedup
# 77359 variants loaded from .bim file.
# 72057 variants pass filters and QC.

## Convert PLINK.1 files to PLINK.2 binary format
plink2 --bfile preprocess/USTA_Diversity_Study.remap.dedup --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --ref-allele 'force' $equCab3_map 7 3 --real-ref-alleles \
        --make-pgen --sort-vars \
        --output-chr 'chrM' --out preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2

##########################################
## Data exploration
##########################################
mkdir -p inspect
## --check-sex compares sex assignments in the input dataset with those imputed from chrX inbreeding coefficients 
## Preliminary run of --check-sex without removing PAR regions
plink2 --pfile preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
      --check-sex max-female-xf=0.2 min-male-xf=0.8 \
      --output-chr 'chrM' --out inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex
tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck | tr ' ' '\t' | cut -f3-5 | sort | uniq -c
#    285 1       1       OK
#    248 2       2       OK
#     43 2       NA      PROBLEM

## Histogram of X chromosome inbreeding coefficients (output of preliminary --check-sex)
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck) > inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck.histo

## Define Par regions and remove them (using high confidence males)
awk '$6 > 0.95 {print $1, $2}' inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.sex.sexcheck > inspect/hiConf_males.txt
#plink2 --pfile preprocess/USTA_Diversity_Study --keep inspect/hiConf_males.txt --het --out inspect/male_het
cat preprocess/USTA_Diversity_Study.sex  | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"0"}' > preprocess/USTA_Diversity_Study.NoSex 
plink2 --pfile preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
       --keep inspect/hiConf_males.txt --update-sex preprocess/USTA_Diversity_Study.NoSex --geno-counts --out inspect/freqx_male
head -n1  inspect/freqx_male.gcount >  inspect/freqx_male.X.gcount
grep "^X" inspect/freqx_male.gcount >>  inspect/freqx_male.X.gcount
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($1=="chrX")a[$3]=$2;next}{$1=a[$2];print $0;}' preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.pvar inspect/freqx_male.X.gcount > inspect/freqx_male.X_ann.gcount
## Heterozygosity is seen until 2063653 which match our expectations (The PAB location on the X chromosome of EquCab2 is located at 1,175,430 bp)
## There is no detected coordinates for the tail PAR, thus I will use the position of the last marker + 1

## --check-sex with removal of PAR and noisy regions
awk 'BEGIN{FS="\t"}{if($1<2063653 || $6>5)print $2}' inspect/freqx_male.X_ann.gcount >  inspect/PAR_and_noise.list
plink2 --pfile preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
    --exclude inspect/PAR_and_noise.list \
    --check-sex max-female-xf=0.2 min-male-xf=0.8 \
    --output-chr 'chrM' --out inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex

tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck | tr ' ' '\t' | cut -f3-5 | sort | uniq -c
#    285 1       1       OK
#    247 2       2       OK
#     44 2       NA      PROBLEM
##  285 (male) and 247 (female) and 44 (NA). (Females are lost due to higher inbreeding)

## Histogram of X chromosome inbreeding coefficients (output of final --check-sex)
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck) > inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck.histo


## --het, --missing, --freq, --hardy
plink2 --pfile preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
      --het --missing --freq --hardy 'midp'  \
      --output-chr 'chrM' --out inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore

#--freq: Allele frequencies (founders only) written to inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.afreq .
#--missing: Sample missing data report written to inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.smiss .
#--missing: Variant missing data report written to inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.vmiss .
#--hardy midp: Autosomal Hardy-Weinberg report (founders only) written to inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy .
#--hardy midp: chrX Hardy-Weinberg report (founders only) written to inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy.x .
#Excluding 3581 variants on non-autosomes from --het.
#--het: done.
#Warning: 3937 variants skipped because they were monomorphic. 
#use --read-freq to provide more accurate allele frequency estimates.
#--het: Results written to inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het .

## heterozygosity (Postive estimates indicate high homozygosity while negative estimates indicate low homozygosity.)
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                END { for(i=bmin;i<=bmax;++i){if(i==0) print -1*size,size,a[i]/1;else if(i<0) print (i-1)*size,i*size,a[i]/1;else print i*size,(i+1)*size,a[i]/1 }}'  <(tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het) > inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het.histo
```
-0.1    -0.08   3
-0.08   -0.06   6
-0.06   -0.04   16
-0.04   -0.02   27
-0.02   0.02    114
0.02    0.04    85
0.04    0.06    96
0.06    0.08    84
0.08    0.1     60
0.1     0.12    30
0.12    0.14    22
0.14    0.16    12
0.16    0.18    11
0.18    0.2     6
0.2     0.22    2
0.22    0.24    1
0.24    0.26    1
```
paste inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het | cut -f1,2,4,6,9-12 | less
paste inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2_noPAR.sex.sexcheck inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het | cut -f1,2,4,6,9-12 | awk -F"\t" '{if($3=="NA")print}' | sort -t $'\t' -k4,4g
#cat inspect/AxiomGT1v2.explore.het | awk '{if($6>0.3)print}' >  inspect/AxiomGT1v2.explore.het.highHomo ## high homozygosity (i.e. low heterozygosity)
#cat inspect/AxiomGT1v2.explore.het | awk '{if(NR==1)print}{if($6<-0.3)print}' >  inspect/AxiomGT1v2.explore.het.lowHomo ## low homozygosity


## HWE
awk 'BEGIN{OFS="\t";}{ if($10<1e-50)a["1e-50 or less"]++;
                       else if($10<1e-40)a["1e-40:1e-50"]++; else if($10<1e-30)a["1e-30:1e-40"]++; \
                       else if($10<1e-20)a["1e-20:1e-30"]++; else if($10<1e-10)a["1e-10:1e-20"]++; \
                       else if($10<1e-9)a["1e-9:1e-10"]++; else if($10<1e-8)a["1e-8:1e-9"]++; \
                       else if($10<1e-7)a["1e-7:1e-8"]++; else if($10<1e-6)a["1e-6:1e-7"]++; \
                       else if($10<1e-5)a["1e-5:1e-6"]++; else if($10<1e-4)a["1e-4:1e-5"]++; \
                       else if($10<0.001)a["1e-3:1e-4"]++; else if($10<0.01)a["1e-2:1e-3"]++; \
                       else a["0.01 or more"]++; } \
                 END { for(i in a) print i,a[i] }'  <(tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy) | sort -g > inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy.histo

```
1e-50 or less   17
1e-40:1e-50     11
1e-30:1e-40     32
1e-20:1e-30     81
1e-10:1e-20     604
1e-9:1e-10      197
1e-8:1e-9       264
1e-7:1e-8       422
1e-6:1e-7       582
1e-5:1e-6       718
1e-4:1e-5       1253
1e-3:1e-4       2004
1e-2:1e-3       4043
0.01 or more    58274
```

## Use this to further explore variants with extreme deviation from HWE:
cat inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy | awk '{if(NR==1)print}{if($10<1e-50)print}' >  inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.hardy.lowHWE



## Variants with very low MAF
# Here are 2 different resolutions for a histogram of MAF:
awk -v size=0.01 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.afreq) > inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.frq.histo
awk -v size=0.001 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.afreq) > inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.frq.histo2


## Final filtering based on missingness, MAF, and HWE
mkdir -p filtered
plink2 --pfile preprocess/USTA_Diversity_Study.remap.refAlleles.dedup.plink2 \
      --hwe 1e-6 'midp' --geno 0.05 --mind 0.05 --maf 0.01 \
      --real-ref-alleles --make-pgen --output-chr 'chrM' --out filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered

#0 samples removed due to missing genotype data (--mind).
#--geno: 493 variants removed due to missing genotype data.
#--hwe midp: 2144 variants removed due to Hardy-Weinberg exact test (founders only).
# 8178 variants removed due to allele frequency threshold(s)
# 61242 variants remaining after main filters.

## Check final genotyping rate
plink2 --pfile filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered --genotyping-rate \
       --out filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered.genotyping_rate ## Total (hardcall) genotyping rate is 0.997762.

## Convert back to PLINK1 binary format for compatibility with other tools
plink2 --pfile filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered \
      --real-ref-alleles --make-bed --output-chr 'chrM' --out filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered


## LD pruning to get independent variants for diversity calculations (& and output as PLINK1 binary format)
mkdir -p LD_pruned
plink2 --pfile filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered \
       --indep-pairwise 100kb 0.8 \
       --real-ref-alleles --output-chr 'chrM' --out LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered.LD_lst ## 12798/61242 variants removed

plink2 --pfile filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered \
       --extract LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered.LD_lst.prune.in \
       --real-ref-alleles --make-bed --output-chr 'chrM' --out LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune ## 48444 variants remaining


## Explore the LD-pruned dataset
plink2 --bfile LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --het --missing --freq --hardy 'midp'  \
      --output-chr 'chrM' --out inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune.explore
paste inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.explore.het inspect/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune.explore.het | less ## F (i.e., measurement of inbreeding) decrease after pruning

## Check final genotyping rate of the LD-pruned dataset
plink2 --bfile LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --genotyping-rate --out LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune.genotyping_rate ## Total (hardcall) genotyping rate is 0.997687.

## Convert to VCF format
plink2 --pfile filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered \
      --real-ref-alleles --export vcf id-paste=iid --output-chr 'chrM' --out filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered

plink2 --pfile filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered \
       --extract LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink2.filtered.LD_lst.prune.in \
       --real-ref-alleles --export vcf id-paste=iid --output-chr 'chrM' --out LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.LD_prune 

############## Stats on diversity ##################
mkdir -p divStats
pl1_filtered="filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered"
pl1_pruned="LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune"

vcf_filtered="filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.vcf"
vcf_pruned="LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.LD_prune.vcf"


# check ref alleles and positions
ref="../Horse_parentage_SNPs/equCab3/download/equCab3.fa"

grep -v '^##chrSet' $vcf_filtered | grep -E "^#|^chr" | bgzip --output $vcf_filtered.test.gz
tabix $vcf_filtered.test.gz
bcftools norm -c ws -f $ref $vcf_filtered.test.gz 1> $vcf_filtered.test.check.vcf 2> $vcf_filtered.test.check.log
## Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      61178/0/0/0/0/0/0
## REF/ALT total/modified/added:   61178/0/0

grep -v '^##chrSet' $vcf_pruned | grep -E "^#|^chr" | bgzip --output $vcf_pruned.test.gz
tabix $vcf_pruned.test.gz
bcftools norm -c ws -f $ref $vcf_pruned.test.gz 1> $vcf_fvcf_prunediltered.test.check.vcf 2> $vcf_pruned.test.check.log
## Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      48382/0/0/0/0/0/0
## REF/ALT total/modified/added:   48382/0/0



## Clean chr names for compatibility (This dataset is not used for now)
#plink --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr '0' \
#        --chr 1-31, x \
#        --make-bed \
#        --output-chr 'chrM' --out "$pl1_pruned".clean

##########################################
## 1. Expected and observed heterozygosity and inbreeding coefficient
##########################################
## An inbreeding coefficient (COI) is a measure of the probability that an individual will have two copies of an allele that are identical by descent from a common ancestor. 
## A higher COI means more predictability of traits but also a greater risk of genetic health problems due to inbreeding depression
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --het 'cols=fid,hom,het,nobs,f' \
    --out divStats/filtered.LD_prune.het_stats
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                END { for(i=bmin;i<=bmax;++i){if(i==0) print -1*size,size,a[i]/1;else if(i<0) print (i-1)*size,i*size,a[i]/1;else print i*size,(i+1)*size,a[i]/1 }}'  <(tail -n+2 divStats/filtered.LD_prune.het_stats.het) > divStats/filtered.LD_prune.het_stats.het.histo

rclone -v copy divStats --include "filtered.LD_prune.het_stats.het*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/het_and_COI/" --drive-shared-with-me


#plink --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
#    --het --ibc \
#    --out divStats/filtered.LD_prune.ibc_stats

##########################################
## 2. PCA
##########################################
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --autosome --pca \
       --output-chr 'chrM' --out divStats/filtered.LD_prune.pca

pca_prefix="divStats/filtered.LD_prune.pca"
rclone -v copy divStats --include "filtered.LD_prune.pca.eigen*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

Rscript -e 'args=(commandArgs(TRUE));'\
'val <- read.table(paste(args[1],"eigenval",sep="."));'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = args[2]);'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "principal Component", ylab = "Variance explained in %");'\
'dev.off();' "$pca_prefix" "divStats/Var_PCs.jpg"
rclone -v copy divStats/Var_PCs.jpg "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t";a["IID"]="sex"}NR==FNR{if($5==1)a[$2]="male";else a[$2]="female";next}{print $0,a[$2]}' $pl1_pruned.fam $pca_prefix.eigenvec > $pca_prefix.eigenvec.wSex
eigenvec_suffix="wSex"; color_column="sex"; out_png="divStats/pca_plot_sex.png";
Rscript scripts/pca_6plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

cat $docs/USTA_Gait_BookSize_Assignments_Sex_Added.csv | tr ' ' '_' | tr ',' '\t' > $docs/USTA_Gait_BookSize_Assignments_Sex_Added.tsv
awk 'BEGIN{FS=OFS="\t";a["IID"]="Gait"}NR==FNR{a[$1]=$3;next}{print $0,a[$2]}' $docs/USTA_Gait_BookSize_Assignments_Sex_Added.tsv $pca_prefix.eigenvec > $pca_prefix.eigenvec.wGait
eigenvec_suffix="wGait"; color_column="Gait"; out_png="divStats/pca_plot_Gait.png";
Rscript scripts/pca_6plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$1]=$5;next}{print $0,a[$2]}' $docs/USTA_Gait_BookSize_Assignments_Sex_Added.tsv $pca_prefix.eigenvec > $pca_prefix.eigenvec.wBook_Size 
eigenvec_suffix="wBook_Size"; color_column="Book_Size"; out_png="divStats/pca_plot_BookSize.png";
Rscript scripts/pca_6plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 divStats/filtered.LD_prune.het_stats.het) $pca_prefix.eigenvec > $pca_prefix.eigenvec.wCOI
eigenvec_suffix="wCOI"; color_column="COI"; out_png="divStats/pca_plot_inbreeding.png";
Rscript scripts/pca_6plots_scaleColor.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me


## Identify Trotter samples segregating on PC2
cat $pca_prefix.eigenvec | awk 'BEGIN{FS=OFS="\t"}{if($4>0.1)print $2}' | grep -Fwf - $docs/USTA_Gait_BookSize_Assignments_Sex_Added.csv > divStats/Trotters_segregating_on_PC2.csv
rclone -v copy divStats/Trotters_segregating_on_PC2.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
## Identify Pacer samples co-segregating with Trotters on PC1
cat $pca_prefix.eigenvec | awk 'BEGIN{FS=OFS="\t"}{if($3<0)print $2}' | grep -Fwf - $docs/USTA_Gait_BookSize_Assignments_Sex_Added.csv | grep "Pacer" > divStats/Pacers_cosegregating_withTrotters_on_PC1.csv
rclone -v copy divStats/Pacers_cosegregating_withTrotters_on_PC1.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me
## Identify Trotter samples co-segregating with Pacers on PC1
cat $pca_prefix.eigenvec | awk 'BEGIN{FS=OFS="\t"}{if($3>0)print $2}' | grep -Fwf - $docs/USTA_Gait_BookSize_Assignments_Sex_Added.csv | grep "Trotter" > divStats/Trotters_cosegregating_withPacers_on_PC1.csv
rclone -v copy divStats/Trotters_cosegregating_withPacers_on_PC1.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me


## 2A. PCA (Trotter only)
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --keep <(grep "Trotter" preprocess/USTA_Diversity_Study.gait) \
       --autosome --pca \
       --output-chr 'chrM' --out divStats/filtered.LD_prune.Trotter.pca

pca_prefix="divStats/filtered.LD_prune.Trotter.pca"
rclone -v copy divStats --include "filtered.LD_prune.Trotter.pca.eigen*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

Rscript -e 'args=(commandArgs(TRUE));'\
'val <- read.table(paste(args[1],"eigenval",sep="."));'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = args[2]);'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "principal Component", ylab = "Variance explained in %");'\
'dev.off();' "$pca_prefix" "divStats/Var_PCs.Trotter.jpg"
rclone -v copy divStats/Var_PCs.Trotter.jpg "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$1]=$5;next}{print $0,a[$2]}' $docs/USTA_Gait_BookSize_Assignments_Sex_Added.tsv $pca_prefix.eigenvec > $pca_prefix.eigenvec.wBook_Size 
eigenvec_suffix="wBook_Size"; color_column="Book_Size"; out_png="divStats/pca_plot_BookSize.Trotter.png";
Rscript scripts/pca_3plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 divStats/filtered.LD_prune.het_stats.het) $pca_prefix.eigenvec > $pca_prefix.eigenvec.wCOI
eigenvec_suffix="wCOI"; color_column="COI"; out_png="divStats/pca_plot_inbreeding.Trotter.png";
Rscript scripts/pca_3plots_scaleColor.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me


## 2B. PCA (Pacer only)
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
       --keep <(grep "Pacer" preprocess/USTA_Diversity_Study.gait) \
       --autosome --pca \
       --output-chr 'chrM' --out divStats/filtered.LD_prune.Pacer.pca

pca_prefix="divStats/filtered.LD_prune.Pacer.pca"
rclone -v copy divStats --include "filtered.LD_prune.Pacer.pca.eigen*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

Rscript -e 'args=(commandArgs(TRUE));'\
'val <- read.table(paste(args[1],"eigenval",sep="."));'\
'val$varPerc <- val$V1/sum(val$V1);'\
'jpeg(file = args[2]);'\
'plot( x = seq(1:length(val$varPerc)), y = val$varPerc, type = "o",xlab = "principal Component", ylab = "Variance explained in %");'\
'dev.off();' "$pca_prefix" "divStats/Var_PCs.Pacer.jpg"
rclone -v copy divStats/Var_PCs.Pacer.jpg "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t";a["IID"]="Book_Size"}NR==FNR{a[$1]=$5;next}{print $0,a[$2]}' $docs/USTA_Gait_BookSize_Assignments_Sex_Added.tsv $pca_prefix.eigenvec > $pca_prefix.eigenvec.wBook_Size 
eigenvec_suffix="wBook_Size"; color_column="Book_Size"; out_png="divStats/pca_plot_BookSize.Pacer.png";
Rscript scripts/pca_3plots.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

awk 'BEGIN{FS=OFS="\t";a["IID"]="COI"}NR==FNR{a[$2]=$8;next}{print $0,a[$2]}' <(tail -n+2 divStats/filtered.LD_prune.het_stats.het) $pca_prefix.eigenvec > $pca_prefix.eigenvec.wCOI
eigenvec_suffix="wCOI"; color_column="COI"; out_png="divStats/pca_plot_inbreeding.Pacer.png";
Rscript scripts/pca_3plots_scaleColor.R "$pca_prefix" "$eigenvec_suffix" "$color_column" "$out_png"
rclone -v copy "$out_png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/PCA/" --drive-shared-with-me

##########################################
## 3. Fst between subpopulations (genders, gait types, and book sizes)
##########################################
## The fixation index can range from 0 to 1, where 0 means complete sharing of genetic material and 1 means no sharing. 
## For values equal to 1(meaning no sharing), scientists say that the populations are fixed.
## Effects of marker type and filtering criteria on QST-FST comparisons: https://pmc.ncbi.nlm.nih.gov/articles/PMC6894560/
for group in sex gait bookSize; do
    plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --pheno preprocess/USTA_Diversity_Study.$group \
        --fst 'PHENO1' 'blocksize=2000' \
        --output-chr 'chrM' --out divStats/filtered.LD_prune.fst_$group
done

find divStats/filtered.LD_prune.fst_*.summary -maxdepth 1 -type f | grep -v "\.x\." | xargs cat > divStats/autosomal.fst.summary
rclone -v copy divStats/autosomal.fst.summary "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me

#find divStats/filtered.LD_prune.fst_*.summary -maxdepth 1 -type f | grep "\.x\." | xargs cat > divStats/chrX.fst.summary
#rclone -v copy divStats/chrX.fst.summary "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/" --drive-shared-with-me

##########################################
# 4. Runs of homozygosity (ROH)
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

#for group in sex gait bookSize; do
#    plink --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
#        --pheno preprocess/USTA_Diversity_Study.$group \
#        --homozyg group-verbose 'extend' \
#        --output-chr 'chrM' --out divStats/filtered.LD_prune.roh_$group
#done


group="gait"
case=$(head -n1 preprocess/USTA_Diversity_Study.$group | cut -f3)

##########################################
## 4A. ROH using Plink (Pruned dataset)
##########################################
plink --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno preprocess/USTA_Diversity_Study.$group $case \
        --homozyg 'extend' \
        --output-chr 'chrM' --out divStats/filtered.LD_prune.roh_$group


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
    { printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/NR, sum5/NR, sum6/NR }' divStats/filtered.LD_prune.roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 16.64 
##Average of the total length of runs (kb) across all samples: 175,509.79
##Average of the average length of runs (KBAVG) across all samples: 10,444.63


## Rscript that plots the correlation between "KB" and "KBAVG" from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="divStats/filtered.LD_prune.roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="divStats/filtered.LD_prune.roh_$group.hom"
Rscript scripts/correlation_plot_multiway.R $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_pruned/" --drive-shared-with-me
#rclone -v copy $out_prefix.correlation_heatmap.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_pruned/" --drive-shared-with-me

##########################################
## 4B. ROH using Plink (Filtered dataset without LD pruning)
##########################################
plink --bfile "$pl1_filtered" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno preprocess/USTA_Diversity_Study.$group $case \
        --homozyg 'extend' \
        --output-chr 'chrM' --out divStats/filtered.not_pruned.roh_$group

awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/NR, sum5/NR, sum6/NR }' divStats/filtered.not_pruned.roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 33.86
##Average of the total length of runs (kb) across all samples: 352,620.78
##Average of the average length of runs (KBAVG) across all samples: 10,343.87


## Rscript that plots the correlation between  KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="divStats/filtered.not_pruned.roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="divStats/filtered.not_pruned.roh_$group.hom"
Rscript scripts/correlation_plot_multiway.R $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_filtered/" --drive-shared-with-me
#rclone -v copy $out_prefix.correlation_heatmap.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/plink_filtered/" --drive-shared-with-me

##########################################
## 4C. ROH using Plink (Filtered dataset without LD pruning (with group option))
##########################################
plink --bfile "$pl1_filtered" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
        --make-pheno preprocess/USTA_Diversity_Study.$group $case \
        --homozyg group 'extend' \
        --homozyg-window-snp 20 \
        --output-chr 'chrM' --out divStats/filtered.not_pruned.group_roh_$group

awk 'NR > 1{ sum4 += $4; sum5 += $5; sum6 += $6 } END \
    { printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum4/NR, sum5/NR, sum6/NR }' divStats/filtered.not_pruned.group_roh_$group.hom.indiv
##Average Number of runs of homozygosity (NSEG) : 35.40
##Average of the total length of runs (kb) across all samples: 363,484.94
##Average of the average length of runs (KBAVG) across all samples: 10,194.10

## Rscript that plots the correlation between  KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="divStats/filtered.not_pruned.group_roh_$group.hom.indiv" ## to read KB and KBAVG
het_stats="divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="divStats/filtered.not_pruned.group_roh_$group.hom"
Rscript scripts/correlation_plot_multiway.R $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_filtered_gp/" --drive-shared-with-me
#rclone -v copy divStats/$out_prefix.correlation_heatmap.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/plink_filtered_gp/" --drive-shared-with-me

##########################################
# 4D. ROH using bcftools/roh (Filtered dataset without LD pruning)
##########################################
## Prepare VCF for BEAGLE phasing and bcftools roh
grep -v '^##chrSet' $vcf_filtered | grep -E "^#|^chr" | grep -v "^chrX" | bgzip --output $vcf_filtered.auto.gz
tabix $vcf_filtered.auto.gz

## for some reason, we still have 28 duplicate SNPs after previous deduplication steps, so we need to remove them here again
bcftools norm \
  --rm-dup exact \
  -Oz \
  -o $vcf_filtered.norm.vcf.gz \
  $vcf_filtered.auto.gz ## Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      58178/0/0/0/0/28/0
tabix -p vcf $vcf_filtered.norm.vcf.gz

## Run BEAGLE
beagle gt=$vcf_filtered.norm.vcf.gz out=$vcf_filtered.norm.phased nthreads=10
# Effective population size (Ne) is the number of individuals in an idealized population that would experience the same amount of genetic drift or inbreeding as the real, observed population. 
# we should provide this number as an input to Beagle when imputing few samples in the mating app.
grep "Estimated ne" $vcf_filtered.norm.phased.log | awk -F":" '{a+=$2}END{print "Ave. Estimated ne:",a/NR}' # Ave. Estimated ne: 2789.03
tabix -p vcf $vcf_filtered.norm.phased.vcf.gz

## Assess change in genotyping rate
plink2 --vcf $vcf_filtered.norm.vcf.gz --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --genotyping-rate --out $vcf_filtered.norm.genotyping_rate ## Total (hardcall) genotyping rate is 0.997935.
plink2 --vcf $vcf_filtered.norm.phased.vcf.gz --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
      --genotyping-rate --out $vcf_filtered.norm.phased.genotyping_rate ## Total (hardcall) genotyping rate is 1.

## Run bcftools roh
bcftools roh -G30 --estimate-AF - $vcf_filtered.norm.phased.vcf.gz -o divStats/roh_out.txt
##Number of target samples: 576
##Number of --estimate-AF samples: 576
##Number of sites in the buffer/overlap: unlimited
##Number of lines total/processed: 58178/58178
##Number of lines filtered/no AF/no alt/multiallelic/dup: 0/0/0/0/0

## Visualize the raw ROH results 
#roh-viz -i divStats/roh_out.txt -v $vcf_filtered.gz -o divStats/roh_viz.html
#rclone -v copy divStats/roh_viz.html "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/" --drive-shared-with-me

grep -E "^RG|^#" divStats/roh_out.txt > divStats/roh_out_RG.txt
## Summary stats by RG
awk 'BEGIN{print "IID\tNSEG\tKB\tKBAVG"} $1=="RG"{n[$2]++; sum[$2]+=$6} END{for (s in n) printf "%s\t%d\t%.2f\t%.2f\n", s, n[s], sum[s]/1000, (sum[s]/1000)/n[s]}' divStats/roh_out_RG.txt > divStats/roh_summary_by_RG.txt
awk 'NR > 1{ sum2 += $2; sum3 += $3; sum4 += $4 } END \
    { printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum2/NR, sum3/NR, sum4/NR }' divStats/roh_summary_by_RG.txt
##Average Number of runs of homozygosity (NSEG) : 86.58
##Average of the total length of runs (kb) across all samples: 454,816.26
##Average of the average length of runs (KBAVG) across all samples: 5,229.99

## filtration to match the PLINK quality suggestions 
#Minimum ROH length (--homozyg-kb) 1000 kb
awk '/^#/ || $6 >= 1000000' divStats/roh_out_RG.txt > divStats/roh.L1.txt
#Minimum number of SNPs in ROH (--homozyg-snp) 50
awk '/^#/ || $7 >= 50' divStats/roh.L1.txt > divStats/roh.L2.txt
#Quality scores
awk -v size=2 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(grep -v "^#" divStats/roh.L2.txt) > divStats/roh.L2.histo 
awk '/^#/ || $8 >= 20' divStats/roh.L2.txt > divStats/roh.L3.txt

## Visualize the ROH results after QC filtration
#roh-viz -i divStats/roh.L3.txt -v $vcf_filtered.gz -o divStats/roh.L3_viz.html
#rclone -v copy divStats/roh.L3_viz.html "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/" --drive-shared-with-me

## Summary stats by RG after QC filtration
awk 'BEGIN{print "IID\tNSEG\tKB\tKBAVG"} $1=="RG"{n[$2]++; sum[$2]+=$6} END{for (s in n) printf "%s\t%d\t%.2f\t%.2f\n", s, n[s], sum[s]/1000, (sum[s]/1000)/n[s]}' divStats/roh.L3.txt > divStats/roh_summary_by_RG_L3.txt
awk 'NR > 1{ sum2 += $2; sum3 += $3; sum4 += $4 } END \
    { printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum2/NR, sum3/NR, sum4/NR }' divStats/roh_summary_by_RG_L3.txt
##Average Number of runs of homozygosity (NSEG) : 54.60
##Average of the total length of runs (kb) across all samples: 411,787.308
##Average of the average length of runs (KBAVG) across all samples: 7,492.35

## Stratify the file by the gait type
## Pacer
awk 'BEGIN{gait["IID"]="gait"}FNR==NR{gait[$2]=$3;next} {print $0,gait[$1]}' preprocess/USTA_Diversity_Study.gait divStats/roh_summary_by_RG_L3.txt > divStats/roh.L3_gait.txt
awk '{ sum2 += $2; sum3 += $3; sum4 += $4 } END \
    { printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum2/NR, sum3/NR, sum4/NR }' <(grep "Pacer" divStats/roh.L3_gait.txt)
##Average Number of runs of homozygosity (NSEG) : 50.76
##Average of the total length of runs (kb) across all samples: 382,228.83
##Average of the average length of runs (KBAVG) across all samples: 7,515.25

## Trotter
awk '{ sum2 += $2; sum3 += $3; sum4 += $4 } END \
    { printf "Average Number of runs of homozygosity (NSEG) : %.2f\n \
    Average of the total length of runs (kb) across all samples: %.2f\n \
    Average of the average length of runs (KBAVG) across all samples: %.2f\n", \
    sum2/NR, sum3/NR, sum4/NR }' <(grep "Trotter" divStats/roh.L3_gait.txt)
##Average Number of runs of homozygosity (NSEG) : 58.62
##Average of the total length of runs (kb) across all samples: 442,775.16
##Average of the average length of runs (KBAVG) across all samples: 7,495.47

## Rscript that plots the correlation between KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
roh_indiv="divStats/roh_summary_by_RG_L3.txt" ## to read KB and KBAVG
het_stats="divStats/filtered.LD_prune.het_stats.het"        ## to read O(HET), E(HET), and F
out_prefix="divStats/filtered.not_pruned.roh_summary_by_RG_L3"
Rscript scripts/correlation_plot_multiway.R $roh_indiv $het_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

## A per-base consensus ROH where ≥25% of samples are in ROH filtered by minimum size 500 kb and stratified by gait type
## Note: we may consider applying Smoothing or LOESS to the per-base coverage data to reduce noise before identifying consensus ROH regions
#roh_RG=divStats/roh_out_RG
roh_RG=divStats/roh.L3
# 1. Convert RG output → BED format
awk 'BEGIN{OFS="\t"} $1=="RG" {print $3, $4-1, $5, $2}' "${roh_RG}.txt" > "${roh_RG}.bed"
# 2. Ensure ROHs from the same sample do not double-count
cut -f4 "${roh_RG}.bed" | sort -u | while read S; do
  awk -v s="$S" '$4==s' "${roh_RG}.bed" | sort -k1,1 -k2,2n | bedtools merge -i - | awk -v s="$S" 'BEGIN{OFS="\t"}{print $1,$2,$3,s}'
done | sort -k1,1 -k2,2n > "${roh_RG}.merged_per_sample.wholePop.bed"
grep "Trotter" preprocess/USTA_Diversity_Study.gait | cut -f2 | grep -f - "${roh_RG}.merged_per_sample.wholePop.bed" > "${roh_RG}.merged_per_sample.Trotter.bed"
grep "Pacer" preprocess/USTA_Diversity_Study.gait | cut -f2 | grep -f - "${roh_RG}.merged_per_sample.wholePop.bed" > "${roh_RG}.merged_per_sample.Pacer.bed"
# 3. Calculate per-base ROH frequency (i.e., how many samples are in ROH at each base position)
reference_fai=$HOME/Equine80select_remapper/equCab3/equCab3_genome.fa.fai
awk '$1 ~ /^[0-9]+$/' $reference_fai | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2}' > divStats/autosomes.genome
for rg in "wholePop" "Trotter" "Pacer"; do
    bedtools genomecov -i "${roh_RG}.merged_per_sample.${rg}.bed" -g divStats/autosomes.genome -bg > "${roh_RG}.per_base_coverage.${rg}.bed"
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
pct=25
for rg in "wholePop" "Trotter" "Pacer"; do
  num_samples=$(cut -f4 "${roh_RG}.merged_per_sample.${rg}.bed" | sort -u | wc -l)
  threshold=$(echo "$pct * $num_samples / 100" | bc -l)
  ## Find consensus before smoothing
  awk -v threshold=$threshold 'BEGIN{OFS="\t"} $4 >= threshold {print}' "${roh_RG}.per_base_coverage.${rg}.bed" > "${roh_RG}.consensus_${pct}pct.${rg}.bed"
  bedtools merge -i "${roh_RG}.consensus_${pct}pct.${rg}.bed" -c 4 -o mean | awk 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=0.5)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${rg}.bed"

  # Summary stats of consensus ROH regions
  echo "==== consensus ROH in ≥${pct}% of ${rg} samples ======"
  awk -v rg="$rg" 'BEGIN{maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
    {print rg,"\nNo of ROH regions:",NR,\
    "\nTotal length of consensus ROH regions (Mbp):", sumLen,\
    "\nAverage ROH length (Mbp):",sumLen/NR,\
    "\nMax no of samples in consensus:", maxConsen,\
    "\nAverage no of samples in consensus", sumSamples/NR}' "${roh_RG}.consensus_${pct}pct.merged.${rg}.bed"

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
  bedtools merge -i "${roh_RG}.consensus_${pct}pct.${rg}.smoothed.bed" -c 4 -o mean | awk 'BEGIN{FS=OFS="\t"}{size=($3-$2)/1000000;if(size>=0.5)print $0,size}' > "${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed"
  
  # upload bed files
  # pause for now to save space
  #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${rg}.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
  #rclone -v copy ${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me
  
  # Summary stats of consensus ROH regions
  echo "==== consensus smoothed ROH in ≥${pct}% of ${rg} samples ======"
  awk -v rg="$rg" 'BEGIN{maxConsen=0;sumSamples=0;sumLen=0;} {if(maxConsen<$4)maxConsen=$4; sumSamples += $4; sumLen += $5} END \
    {print rg,"\nNo of ROH regions:",NR,\
    "\nTotal length of consensus ROH regions (Mbp):", sumLen,\
    "\nAverage ROH length (Mbp):",sumLen/NR,\
    "\nMax no of samples in consensus:", maxConsen,\
    "\nAverage no of samples in consensus", sumSamples/NR}' "${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed"
done > ${roh_RG}.consensus_${pct}pct.summary.txt
rclone -v copy ${roh_RG}.consensus_${pct}pct.summary.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me


## intersect ROH regions of each sample aganist the consensus wholePop ROH regions
roh_RG="divStats/roh.L3"
bed_perSample="${roh_RG}.merged_per_sample.wholePop.bed"
rg="wholePop"
consensus_pct=${roh_RG}.consensus_${pct}pct.merged.${rg}.smoothed.bed
consensus_size=$(awk 'BEGIN{sum=0} {sum+=($3-$2)} END {print sum}' ${consensus_pct})
echo -e "IID\tTotal_ROH_in_Consensus_region(bp)\tPercent_of_Consensus_ROH" > ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt
cut -f4 "${bed_perSample}" | sort -u | while read S; do
  awk -v s="$S" '$4==s' "${bed_perSample}" | sort -k1,1 -k2,2n | bedtools intersect -a stdin -b "${consensus_pct}" | awk -v s="$S" -v cs="$consensus_size" 'BEGIN{OFS="\t"}{size+=($3-$2)} END {print s, size, (size/cs)*100}'
done >> ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt
rclone -v copy ${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/ROH/bcftools/" --drive-shared-with-me

############################################
## 4E. Studying ROH using Howard et al. (2016) approach
############################################
mkdir -p $HOME/genDiv/rep_ROHRM
rohrm_dir="rep_ROHRM"

mamba create -n GENpy scikit-allel numpy pandas
conda activate GENpy

## ROH Analysis with Sub-populations and Phenotypes
## Having a headerless tab-separated file with 3 columns: The 2nd column has subject ids matching the VCF and the 3rd column has the a binary phenotype, let us do the following:
## 1. ROH Calling (Per Individual): We will implement a scanner that checks hap1 == hap2. Any contiguous stretch of matching haplotypes longer than the cutoff (e.g., 1 Mb) is flagged as an ROH.
## 2. Island Detection: We will calculate the frequency of ROHs at every SNP, find the "Top 5%" cutoff, and merge contiguous high-frequency SNPs into "Islands".
##    i.e., an "ROH Island" is defined strictly as a contiguous block of SNPs where every single SNP is in the Top 5% of frequencies.
## 3. Phenotype Integration: the analysis will be repeated for sub-populations defined in the phenotype file).
## 4. make a plot to show the ROH frequency across the genome for the whole population and each sub-population.
## 5. Calc the average (±SD) proportion of the genome in a ROH for the whole population and each sub-population.

roh_mb_cutoff=1.0  # in Megabases (Mb)
phenotypes="preprocess/USTA_Diversity_Study.$group"
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

############################################
## 5. F_ROH statistic (currently calculated based on bacftools roh)
############################################
## Standard practice is to calculate F_{ROH} statistics on all valid ROHs (>1Mb), while restricting "Islands" (signatures of selection) to only the most robust regions.
## per-sample F_ROH = (sum length of ROH for that individual) / (total autosomal genome length).

## calculate the effective autosomal genome length
awk '{print $1"\t"$4}' "$pl1_filtered".bim | grep "^chr" | grep -v "^chrX" > "$pl1_filtered".snp_pos.txt
aut_len=$(sort -k1,1 -k2,2n "$pl1_filtered".snp_pos.txt | \
        awk '{if ($1 == prev_chr) { gap = $2 - prev_pos; \
              if(gap > 0) {if (gap > 1000000) gap = 1000000; total += gap; }}\
              prev_chr=$1; prev_pos=$2} END {print total}') ## 2,261,547,402

awk -v aut_len=$aut_len 'BEGIN{FS=OFS="\t";}NR==1{print $0,"F_ROH";next} {print $0, ($3*1000)/aut_len}' divStats/roh_summary_by_RG_L3.txt > divStats/roh_summary_by_RG_L3_Froh.txt
awk -v size=0.02 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 divStats/roh_summary_by_RG_L3_Froh.txt) > divStats/roh_summary_by_RG_L3_Froh.histo 
rclone -v copy divStats  --drive-shared-with-me --include "roh_summary_by_RG_L3_Froh.*" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/"

awk '{if($5>0.3)print $0}' divStats/roh_summary_by_RG_L3_Froh.txt | tr '\t' ',' > divStats/roh_high.csv
rclone -v copy divStats/roh_high.csv  --drive-shared-with-me "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/"



############################################
## 6. Relatedness work
############################################


############################################
## Genomic Relatedness using Plink
############################################
## Create a standard GRM matrix 
conda activate grGWAS
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --pheno $phenotypes \
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
```
Feature,            Standard GRM (VanRaden),                    ROH GRM (Howard et al.)
Input Data,         "Unphased Genotypes (0, 1, 2)",             Phased Haplotypes (0|0, 0|1, 1|0, 1|1)
Resolution,         Single Nucleotide (SNP),                    Multi-Megabase Window (Window)
Matching Logic,     Allele Sharing (IBS),                       Exact String Matching (IBD)
Sensitivity,        Very tolerant of mutation/recombination.,   Very strict. One mismatch breaks the link.
Biological Signal,  Captures Deep/Ancient Relatedness.,         Captures Recent Relatedness.
Diagonal,           Heterozygosity-based Inbreeding.,           ROH-based Inbreeding (FROH​).
```
roh_mb_cutoff=1.0  # in Megabases (Mb)
roh_threshold=3.0  # in Standard Deviations (SD)
conda activate GENpy
python $scripts/ROHRM_Creator.py $vcf_filtered.norm.phased.vcf.gz $roh_mb_cutoff $roh_threshold "$rohrm_dir" > $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold.log
## Stats: Mean=27.91, SD=5.11
## Cutoff Calculated: 13 SNPs
## Windows filtered: 57309 -> 56971 remaining.

############################################
## Compare Genomic Relatedness and ROH-based Relatedness 
############################################
## compares the ROH-based matrix against a Standard vanraden Genomic Relationship Matrix ($G_{STD}$) 
## to demonstrate that $G_{ROH}$ captures different genetic signals (recent ancestry vs. deep ancestry).

## Old version comparing with vanraden GRM implementation. Also, no QC of ID matching and simplified plotting
## You may ignore this and use the analysis_comparison.py script below
#python $scripts/analysis_comparison_vanradenGRM.py $vcf_filtered.norm.phased.vcf.gz $rohrm_dir/ROHRM.rohMinSize_1.0.rohThreshold_3.0.txt "$rohrm_dir"
##Correlation between Inbreeding estimates: 0.3963
##Correlation between Relationship estimates: 0.9286
##Plot saved to 'ROHRM_vs_VanradenGRM_Comparison.png'
#rclone -v copy $rohrm_dir/ROHRM_vs_VanradenGRM_Comparison.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## New version comparing with plink2 GRM implementation. Also, includes QC of ID matching and improved plotting
## python analysis_comparison.py <ROH_Prefix> <Std_Prefix> <Phenotypes_File> <output_dir>
python $scripts/analysis_comparison.py $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold $rohrm_dir/filtered.LD_prune.GRM_$group $phenotypes "$rohrm_dir"
#Correlation (Inbreeding): r = -0.1686
#Correlation (Relationships): r = 0.9288
#Correlations per Subpopulation:
#  > Trotter: r = -0.1595
#  > Pacer: r = 0.1842

# center column 7 (Difference) around the column's mean and save to new column centered_Kinship_diff
awk -F, 'BEGIN{OFS=FS=","} NR==1{hdr=$0; next} {sum+=$7; n++; lines[n]=$0; vals[n]=$7} \
         END{ mean = (n?sum/n:0); print hdr, "centered_Kinship_diff"; \
              for(i=1;i<=n;i++) printf "%s%s%.8f\n", lines[i], OFS, vals[i]-mean }' $rohrm_dir/Pairwise_Differences.csv > $rohrm_dir/Pairwise_Differences.csv.tmp && mv $rohrm_dir/Pairwise_Differences.csv.tmp $rohrm_dir/Pairwise_Differences.csv

# generate histograms
awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_Std.histo
awk -v size=0.05 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_ROH.histo
awk -v size=0.01 'BEGIN{FS=",";OFS="\t";bmin=bmax=0}{ b=int($8/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 $rohrm_dir/Pairwise_Differences.csv) > $rohrm_dir/Pairwise_Differences.Kinship_diff.histo

## animal pairs with high positive kinship difference (ROH-based kinship is higher relative to standard kinship)
cat $rohrm_dir/Pairwise_Differences.csv | awk 'BEGIN{FS=","} NR==1{print;next}{if($8>0.1) print}'| head #> high_positive_kinship_diff.csv
## animal pairs with high negative kinship difference (ROH-based kinship is lower relative to standard kinship)
cat $rohrm_dir/Pairwise_Differences.csv | awk 'BEGIN{FS=","} NR==1{print;next}{if($8<-0.15) print}'| head #> high_negative_kinship_diff.csv

cat $rohrm_dir/Pairwise_Differences.csv | grep "Trotter" | grep "Pacer" | awk 'BEGIN{FS=","} {sum+=$8; n++} END{print sum/n}'  # average centered kinship difference between Trotter and Pacer: 0.00185561
cat $rohrm_dir/Pairwise_Differences.csv | grep "Trotter" | grep -v "Pacer" | awk 'BEGIN{FS=","} {sum+=$8; n++} END{print sum/n}'  # average centered kinship difference between Trotters: 0.0222411 (i.e., ROH-based kinship is relatively higher than standard kinship on average for Trotters)
cat $rohrm_dir/Pairwise_Differences.csv | grep -v "Trotter" | grep "Pacer" | awk 'BEGIN{FS=","} {sum+=$8; n++} END{print sum/n}'  # average centered kinship difference between Pacers: -0.0259652 (i.e., standard kinship is relatively higher than ROH-based kinship on average for Pacers)

# upload results
rclone -v copy $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_1Mb.Threshold_3SD/" --drive-shared-with-me
rclone -v copy $rohrm_dir/Pairwise_Differences.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_1Mb.Threshold_3SD/" --drive-shared-with-me
rclone -v copy $rohrm_dir/Inbreeding_Comparison.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_1Mb.Threshold_3SD/" --drive-shared-with-me

# move results to a separate folder
mkdir -p $rohrm_dir/roh_1Mb.Threshold_3SD 
mv $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png $rohrm_dir/Inbreeding_Comparison.csv $rohrm_dir/Pairwise_Differences.* $rohrm_dir/roh_1Mb.Threshold_3SD/
############################################
## test other roh size cutoffs
############################################
roh_mb_cutoff=5.0  # in Megabases (Mb)
roh_threshold=3.0  # in Standard Deviations (SD)
python $scripts/ROHRM_Creator.py $vcf_filtered.norm.phased.vcf.gz $roh_mb_cutoff $roh_threshold "$rohrm_dir" > $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold.log
python $scripts/analysis_comparison.py $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold $rohrm_dir/filtered.LD_prune.GRM_$group $phenotypes "$rohrm_dir"
# Global Correlation (Inbreeding): r = -0.0398
# Global Correlation (Relationships): r = 0.9383
#Correlations per Subpopulation:
#  > Trotter: r = -0.0631
#  > Pacer: r = 0.2412
awk -F, 'BEGIN{OFS=FS=","} NR==1{hdr=$0; next} {sum+=$7; n++; lines[n]=$0; vals[n]=$7} \
         END{ mean = (n?sum/n:0); print hdr, "centered_Kinship_diff"; \
              for(i=1;i<=n;i++) printf "%s%s%.8f\n", lines[i], OFS, vals[i]-mean }' $rohrm_dir/Pairwise_Differences.csv > $rohrm_dir/Pairwise_Differences.csv.tmp && mv $rohrm_dir/Pairwise_Differences.csv.tmp $rohrm_dir/Pairwise_Differences.csv
rclone -v copy $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_5Mb.Threshold_3SD/" --drive-shared-with-me
rclone -v copy $rohrm_dir/Pairwise_Differences.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_5Mb.Threshold_3SD/" --drive-shared-with-me
rclone -v copy $rohrm_dir/Inbreeding_Comparison.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_5Mb.Threshold_3SD/" --drive-shared-with-me
mkdir -p $rohrm_dir/roh_5Mb.Threshold_3SD
mv $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png $rohrm_dir/Inbreeding_Comparison.csv $rohrm_dir/Pairwise_Differences.* $rohrm_dir/roh_5Mb.Threshold_3SD/

roh_mb_cutoff=10.0  # in Megabases (Mb)
roh_threshold=3.0  # in Standard Deviations (SD)
python $scripts/ROHRM_Creator.py $vcf_filtered.norm.phased.vcf.gz $roh_mb_cutoff $roh_threshold "$rohrm_dir" > $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold.log
python $scripts/analysis_comparison.py $rohrm_dir/ROHRM.rohMinSize_$roh_mb_cutoff.rohThreshold_$roh_threshold $rohrm_dir/filtered.LD_prune.GRM_$group $phenotypes "$rohrm_dir"
# Global Correlation (Inbreeding): r = 0.0613
# Global Correlation (Relationships): r = 0.9131
#Correlations per Subpopulation:
#  > Trotter: r = 0.0180
#  > Pacer: r = 0.2604
awk -F, 'BEGIN{OFS=FS=","} NR==1{hdr=$0; next} {sum+=$7; n++; lines[n]=$0; vals[n]=$7} \
         END{ mean = (n?sum/n:0); print hdr, "centered_Kinship_diff"; \
              for(i=1;i<=n;i++) printf "%s%s%.8f\n", lines[i], OFS, vals[i]-mean }' $rohrm_dir/Pairwise_Differences.csv > $rohrm_dir/Pairwise_Differences.csv.tmp && mv $rohrm_dir/Pairwise_Differences.csv.tmp $rohrm_dir/Pairwise_Differences.csv
rclone -v copy $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_10Mb.Threshold_3SD/" --drive-shared-with-me
rclone -v copy $rohrm_dir/Pairwise_Differences.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_10Mb.Threshold_3SD/" --drive-shared-with-me
rclone -v copy $rohrm_dir/Inbreeding_Comparison.csv "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/roh_10Mb.Threshold_3SD/" --drive-shared-with-me
mkdir -p $rohrm_dir/roh_10Mb.Threshold_3SD
mv $rohrm_dir/Robust_Matrix_Comparison_Enhanced.png $rohrm_dir/Inbreeding_Comparison.csv $rohrm_dir/Pairwise_Differences.* $rohrm_dir/roh_10Mb.Threshold_3SD/


## Temp ###########################################
## Compare Inbreeding Coefficients from standard GRM and ROH-based GRM vs Heterozygosity and COIfficient of Inbreeding (COI)
############################################
## compare with het and coi
## Rscript that plots the correlation between KB and KBAVG from .hom.indiv and the difference O(HET) and E(HET), and F columns from .het
RM_diag="rep_ROHRM/roh_1Mb.Threshold_3SD/Inbreeding_Comparison.csv" ## to read F_Standard (comparable to COI "F" measure in 1) and F_ROH (comparable to F_ROH measured in 5)
het_stats="divStats/filtered.LD_prune.het_stats.het"                ## to read O(HET), E(HET), and F
Froh_stats="divStats/roh_summary_by_RG_L3_Froh.txt"          ## to read F_ROH measured in 5
out_prefix="divStats/coi_Froh_rmDiag_correlation"
Rscript scripts/correlation_plot_multiway_v2.R $RM_diag $het_stats $Froh_stats $out_prefix
rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me


roh_RG="divStats/roh.L3"; rg="wholePop"; conShare=${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt;
out_prefix2="divStats/coi_Froh_rmDiag_conShare_correlation"
Rscript scripts/correlation_plot_multiway_v2e.R $RM_diag $het_stats $Froh_stats $conShare $out_prefix2
rclone -v copy $out_prefix2.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

## Resume relatedness work ############################################
############################################
## KING-robust kinship estimator
############################################
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --pheno preprocess/USTA_Diversity_Study.$group \
    --make-king-table 'counts' 'cols=+ibs1' \
    --output-chr 'chrM' --out divStats/filtered.LD_prune.king_$group

kingkin="divStats/filtered.LD_prune.king_$group.kin0"
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($10/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }'  <(tail -n+2 $kingkin) > ${kingkin%.kin0}.histo
rclone -v copy $kingkin "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
rclone -v copy ${kingkin%.kin0}.histo "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

############################################
## calc IBS
############################################
## IBS1= HET1_HOM2 + HET2_HOM1 
## IBS2= N_SNPs - (HETHET + IBS0 + IBS1)
## IBS = (2*IBS2 + IBS1) / (2*N_SNPs)
#awk 'BEGIN{FS=OFS="\t"}NR==1{print $0,"IBS";next}{print $0,(2*$7+$8+$9)/(2*$5)}' $kingkin > ${kingkin}.withIBS
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0,"IBS";next}{ibs1=$8+$9;ibs2=$5-($6+$7+ibs1);print $0,(2*ibs2+ibs1)/(2*$5)}' $kingkin > ${kingkin}.withIBS

## useless
## Plot the correlation between KING-robust kinship and IBS
Rscript $scripts/plot_correlation.R ${kingkin}.withIBS KINSHIP IBS
rclone -v copy correlation_plot_KINSHIP_vs_IBS.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me

############################################
## PCA-based pairwise Euclidean distance
############################################
for pca_prefix in divStats/filtered.LD_prune{\.,\.Trotter\.,\.Pacer\.}pca ;do 
    awk '
    BEGIN {FS=OFS="\t"}
    $1!="#FID" && NF>3 {
        n++
        fid[n]=$1
        iid[n]=$2
        for(c=3;c<=NF;c++) pc[n,c-2] = $c
        npcs = NF - 2
    }
    END {
        print "FID1","IID1","FID2","IID2","PCA_EUCLIDEAN_DIST","DIST_KINSHIP"
        for(i=1;i<=n;i++)
            for(j=1;j<i;j++) {
                dist2=0
                for(c=1;c<=npcs;c++) {
                    d = pc[i,c] - pc[j,c]
                    dist2 += d*d
                }
                dist = sqrt(dist2)
                printf "%s\t%s\t%s\t%s\t%.8f\t%.8f\n",
                    fid[i], iid[i], fid[j], iid[j], dist, exp(-dist2/2)
            }
    }
    ' "$pca_prefix.eigenvec" > "$pca_prefix.pca_pairwise_euclidean.dist"
done

## Useless
## Merge PCA-based Euclidean distance with KING-robust kinship + IBS
kingkin_wIBS=${kingkin}.withIBS
for pca_prefix in divStats/filtered.LD_prune{\.,\.Trotter\.,\.Pacer\.}pca ;do 
    euclDist="$pca_prefix.pca_pairwise_euclidean.dist"
    out_file="$pca_prefix.pca_pairwise_euclidean.dist.withKIN0"
    awk 'BEGIN {FS=OFS="\t"} FNR == NR {
        if ($1 ~ /^#/) next
        # canonical key for pair
        key = ($1":"$2 < $3":"$4) ?
            $1":"$2"|" $3":"$4 :
            $3":"$4"|" $1":"$2
        # store columns 5+ only
        kin0_extra = ""
        for (i = 5; i <= NF; i++)
            kin0_extra = kin0_extra OFS $i

        kin0_data[key] = substr(kin0_extra, 2)   # remove leading OFS
        next
    }
    FNR == 1 {
        print "FID1","IID1","FID2","IID2",
            "PCA_EUCLIDEAN_DIST","KINSHIP_KING_PCA",
            "NSNP","HETHET","IBS0","HET1_HOM2","HET2_HOM1","KINSHIP_PLINK","IBS"
        next
    }
    {
        # build canonical key
        key = ($1":"$2 < $3":"$4) ?
            $1":"$2"|" $3":"$4 :
            $3":"$4"|" $1":"$2

        # fetch .kin0 info (5+ columns)
        extra = (key in kin0_data ? kin0_data[key] : "NA")
        print $1,$2,$3,$4,$5,$6,extra
    }
    ' "$kingkin_wIBS" "$euclDist" > "$out_file"
done

## useless
## Plot the correlation between PCA-based Euclidean distance and KING-robust kinship
for pca_prefix in divStats/filtered.LD_prune{\.,\.Trotter\.,\.Pacer\.}pca ;do 
    out_file="$pca_prefix.pca_pairwise_euclidean.dist.withKIN0"
    Rscript $scripts/plot_correlation.R "$out_file" PCA_EUCLIDEAN_DIST KINSHIP_PLINK
    mv correlation_plot_PCA_EUCLIDEAN_DIST_vs_KINSHIP_PLINK.png "$pca_prefix.correlation_plot_PCA_EUCLIDEAN_DIST_vs_KINSHIP_PLINK.png"
    rclone -v copy "$pca_prefix.correlation_plot_PCA_EUCLIDEAN_DIST_vs_KINSHIP_PLINK.png" "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
done

############################################
## Plot the correlation among ROH-based relatedness, KING-robust kinship, and PCA-based Euclidean distance
############################################
RMs="rep_ROHRM/roh_1Mb.Threshold_3SD/Pairwise_Differences.csv" ## file1 has Kinship_Std and Kinship_ROH columns
kingkin_wIBS="divStats/filtered.LD_prune.king_$group.kin0.withIBS" ## file2 has KINSHIP and IBS columns
for pop in "wholePop" "Trotter" "Pacer";do
    pca_prefix=$(echo divStats/filtered.LD_prune.$pop.pca | sed 's/wholePop\.//')
    euclDist="$pca_prefix.pca_pairwise_euclidean.dist" ## file3 has PCA_EUCLIDEAN_DIST column
    out_prefix="divStats/$pop.relatedness_correlation"
    Rscript scripts/correlation_plot_multiway_v3.R $RMs $kingkin_wIBS $euclDist $out_prefix
    rclone -v copy $out_prefix.pairplot.png "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Relatedness/" --drive-shared-with-me
done

########################################################
## Effective number of alleles (\(A_{e}\)) represents the number of equally frequent alleles required to achieve the same level of expected heterozygosity (\(H_{e}\)) observed in a population

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
## --freq counts: Allele counts (founders only) written to
## LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.LD_prune.freq_stats.afreq
awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.afreq" > "$pl1_pruned.freq_stats.wholePop.afreq.Ae"

plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --pheno preprocess/USTA_Diversity_Study.$group \
    --loop-cats 'PHENO1' --freq \
    --out "$pl1_pruned.freq_stats"

for pop in Trotter Pacer; do
    awk -v pop=$pop 'BEGIN{FS=OFS="\t"} NR==1{print $0,"A_e";next} {p1=$6; p2=1-p1; Ae=1/(p1*p1 + p2*p2); print $0,Ae}' "$pl1_pruned.freq_stats.$pop.afreq" > "$pl1_pruned.freq_stats.$pop.afreq.Ae"
done

## Calculate Mean \(A_{e}\) and standard deviation per population
for pop in wholePop Trotter Pacer; do
    #awk -v pop=$pop 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; n++} END{ mean_Ae = (n?sum_Ae/n:0); sd_Ae = (n?sqrt((sum_Ae^2/n - mean_Ae^2)):0); print "Mean_Ae_in_"pop, mean_Ae, "SD_Ae_in_"pop, sd_Ae }' "$pl1_pruned.freq_stats.$pop.afreq.Ae"
    awk -v pop=$pop 'BEGIN{FS=OFS="\t"} NR==1{next} {sum_Ae+=$NF; sumsq += $NF * $NF; n++} END \
        { if (n > 0) { mean_Ae = sum_Ae/n; sd_Ae = sqrt((sumsq/n - mean_Ae^2)); print "Mean_Ae_in_"pop, mean_Ae, "SD_Ae_in_"pop, sd_Ae } }' "$pl1_pruned.freq_stats.$pop.afreq.Ae"
done
#Mean_Ae_in_wholePop     1.55957 SD_Ae_in_wholePop       0.322746
#Mean_Ae_in_Trotter      1.5311  SD_Ae_in_Trotter        0.33893
#Mean_Ae_in_Pacer        1.55019 SD_Ae_in_Pacer          0.330677


