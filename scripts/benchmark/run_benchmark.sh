#!/usr/bin/env bash
# Benchmark driver for the refactored pipeline.
#
# Given an existing full-run OUTPUT_DIR (the "source"), this script builds a
# 50-sample subset (default 25 Trotter + 25 Pacer), runs the refactored
# genDiversity_per_group.sh × 3 + genDiversity_aggregate.sh against that
# subset, and validates the outputs with scripts/validate_popRefs.sh.
#
# What's subsetted:
#   - PLINK1 filtered BED (${pl1_filtered})
#   - PLINK1 LD-pruned BED (${pl1_pruned})
#   - Phased VCF (${vcf_filtered}.norm.phased.vcf.gz)
#   - samples.${rg}.txt files (wholePop, Trotter, Pacer)
#
# What's regenerated on the subset:
#   - whole-pop KING kinship + IBS augmentation + related (since these are
#     pair statistics and can't be subsetted from the 560-sample originals)
#
# What's copied verbatim from source:
#   - preprocess/USTA_Diversity_Study.{sex,gait,bookSize,gait_bookSize}
#   - preprocess/sample_groups.tsv
#   - divStats/effective_autosomal_genome_length.txt (reference-based constant)
#   - divStats/autosomes.genome (reference-based constant)
#   - Miscellaneous_documents_standardbred/ (the outlier CSVs consume this)
#
# Usage:
#   bash scripts/benchmark/run_benchmark.sh \
#        --source results_20260419_194842 \
#        --target results_benchmark_25x25 \
#        [--trotter-n 25] [--pacer-n 25] [--seed 42]
#
# After it completes, scripts/validate_popRefs.sh reports pass/fail per
# expected file. Non-zero exit = something's missing or malformed; browse
# the target OUTPUT_DIR's run.log for details.

set -eo pipefail

source_dir=""
target_dir=""
trotter_n=25
pacer_n=25
seed=42

while [[ $# -gt 0 ]]; do
    case "$1" in
        --source)     source_dir="$2"; shift 2 ;;
        --target)     target_dir="$2"; shift 2 ;;
        --trotter-n)  trotter_n="$2"; shift 2 ;;
        --pacer-n)    pacer_n="$2"; shift 2 ;;
        --seed)       seed="$2"; shift 2 ;;
        -h|--help)    sed -n '2,32p' "$0"; exit 0 ;;
        *)            echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done
[[ -z "$source_dir" || -z "$target_dir" ]] && { echo "usage: $0 --source SRC --target TGT [--trotter-n 25 --pacer-n 25 --seed 42]" >&2; exit 2; }
[[ ! -d "$source_dir" ]] && { echo "ERROR: source dir missing: $source_dir" >&2; exit 2; }
[[ -e "$target_dir" ]] && { echo "ERROR: target dir already exists: $target_dir (remove or choose another)" >&2; exit 2; }

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
echo "[$(date +%H:%M:%S)] Benchmark: source=$source_dir target=$target_dir Trotter=$trotter_n Pacer=$pacer_n seed=$seed"

# Canonical path stems inside either OUTPUT_DIR
pl1_filtered_stem="filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered"
vcf_filtered_stem="filtered/USTA_Diversity_Study.remap.refAlleles.dedup.vcf.filtered.vcf"
pl1_pruned_stem="LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune"

##############################################################################
# 1. Pick the subset samples (deterministic via --seed)
##############################################################################
mkdir -p "$target_dir/preprocess" "$target_dir/filtered" "$target_dir/LD_pruned" \
         "$target_dir/divStats" "$target_dir/rep_ROHRM"

pick_n() {
    local src="$1" n="$2"
    awk 'NF>=2 {print $1"\t"$2}' "$src" \
        | awk -v s="$seed" 'BEGIN{srand(s); OFS="\t"} {print rand(), $0}' \
        | sort -k1,1n \
        | head -n "$n" \
        | cut -f2-
}

trotter_samples="${source_dir}/preprocess/samples.Trotter.txt"
pacer_samples="${source_dir}/preprocess/samples.Pacer.txt"
[[ -f "$trotter_samples" ]] || { echo "ERROR: missing $trotter_samples" >&2; exit 2; }
[[ -f "$pacer_samples" ]]   || { echo "ERROR: missing $pacer_samples" >&2; exit 2; }

pick_n "$trotter_samples" "$trotter_n" > "$target_dir/preprocess/samples.Trotter.txt"
pick_n "$pacer_samples"   "$pacer_n"   > "$target_dir/preprocess/samples.Pacer.txt"
cat "$target_dir/preprocess/samples.Trotter.txt" "$target_dir/preprocess/samples.Pacer.txt" \
    > "$target_dir/preprocess/samples.wholePop.txt"

echo "[$(date +%H:%M:%S)]   picked $(wc -l < "$target_dir/preprocess/samples.Trotter.txt") Trotter + $(wc -l < "$target_dir/preprocess/samples.Pacer.txt") Pacer = $(wc -l < "$target_dir/preprocess/samples.wholePop.txt") total"

##############################################################################
# 2. Copy reference-based constants and per-sample auxiliary files
##############################################################################
for f in USTA_Diversity_Study.sex USTA_Diversity_Study.gait USTA_Diversity_Study.bookSize USTA_Diversity_Study.gait_bookSize sample_groups.tsv; do
    if [[ -f "${source_dir}/preprocess/${f}" ]]; then
        cp "${source_dir}/preprocess/${f}" "${target_dir}/preprocess/${f}"
    fi
done
cp "${source_dir}/divStats/effective_autosomal_genome_length.txt" "${target_dir}/divStats/"
cp "${source_dir}/divStats/autosomes.genome"                      "${target_dir}/divStats/"
if [[ -d "${source_dir}/Miscellaneous_documents_standardbred" ]]; then
    cp -r "${source_dir}/Miscellaneous_documents_standardbred" "${target_dir}/"
fi

# sample_groups.tsv — regenerated for the 50-sample subset. Matches the
# format shared.sh writes: header + "IID\tgroup" rows, one per sample.
# Samples without a gait label fall back to group=wholePop.
if [[ ! -f "${target_dir}/preprocess/sample_groups.tsv" ]]; then
    {
        echo "# sample_groups.tsv — primary group per reference sample."
        echo "# One row per sample; wholePop membership is implicit."
        echo "# group ∈ {Trotter, Pacer, wholePop}. Samples with no gait label get group=wholePop."
        printf "IID\tgroup\n"
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{gait[$2]=$3;next} {g=gait[$2]; if(g==""||g=="undefined") g="wholePop"; print $2,g}' \
            "${target_dir}/preprocess/USTA_Diversity_Study.gait" "${target_dir}/preprocess/samples.wholePop.txt"
    } > "${target_dir}/preprocess/sample_groups.tsv"
fi

##############################################################################
# 3. Subset PLINK1 filtered, PLINK1 pruned, and the phased VCF
##############################################################################
keep_file="${target_dir}/preprocess/samples.wholePop.txt"

plink2 --bfile "${source_dir}/${pl1_filtered_stem}" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$keep_file" \
    --make-bed --output-chr 'chrM' --out "${target_dir}/${pl1_filtered_stem}"

plink2 --bfile "${source_dir}/${pl1_pruned_stem}" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --keep "$keep_file" \
    --make-bed --output-chr 'chrM' --out "${target_dir}/${pl1_pruned_stem}"

vcf_in="${source_dir}/${vcf_filtered_stem}.norm.phased.vcf.gz"
vcf_out="${target_dir}/${vcf_filtered_stem}.norm.phased.vcf.gz"
bcftools view -S <(cut -f2 "$keep_file") --force-samples "$vcf_in" -Oz -o "$vcf_out"
bcftools index -t "$vcf_out"

##############################################################################
# 4. Regenerate whole-pop KING + IBS + related on the subset
##############################################################################
pl1_pruned="${target_dir}/${pl1_pruned_stem}"
group="gait"
plink2 --bfile "$pl1_pruned" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
    --make-king-table 'counts' 'cols=+ibs1' \
    --output-chr 'chrM' --out "${target_dir}/divStats/filtered.LD_prune.king_${group}"
kingkin="${target_dir}/divStats/filtered.LD_prune.king_${group}.kin0"
awk -v size=0.05 'BEGIN{OFS="\t";bmin=bmax=0}{ b=int($10/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } \
                    END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/1 }' <(tail -n+2 "$kingkin") > "${kingkin%.kin0}.histo"
head -n 1 "$kingkin" > "${target_dir}/divStats/related"
tail -n +2 "$kingkin" | sort -grk10,10 | awk '{if($10>0.177)print}' >> "${target_dir}/divStats/related"
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0,"IBS";next}{ibs1=$8+$9;ibs2=$5-($6+$7+ibs1);print $0,(2*ibs2+ibs1)/(2*$5)}' "$kingkin" > "${kingkin}.withIBS"

##############################################################################
# 5. Run per_group.sh × 3 then aggregate.sh against the target OUTPUT_DIR
##############################################################################
cd "$repo_root"
export OUTPUT_DIR="$target_dir"
for rg in wholePop Trotter Pacer; do
    echo "[$(date +%H:%M:%S)] === per_group.sh $rg ==="
    bash ./genDiversity_per_group.sh "$rg"
done

echo "[$(date +%H:%M:%S)] === aggregate.sh ==="
bash ./genDiversity_aggregate.sh

##############################################################################
# 6. Validate
##############################################################################
echo "[$(date +%H:%M:%S)] === validate_popRefs.sh ==="
bash scripts/validate_popRefs.sh --mode upstream --root "$target_dir"

echo "[$(date +%H:%M:%S)] Benchmark complete: $target_dir"
