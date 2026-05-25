#!/usr/bin/env bash
# Upload-only re-runner: pushes the artifacts already on disk under $OUTPUT_DIR
# to Google Drive, without re-running plink / bcftools / R / Python.
#
# Mirrors every `rclone -v copy` and `upload` call in
# genDiversity_per_group.sh (looped over wholePop / Trotter / Pacer and over
# $ROH_CUTOFFS) and genDiversity_aggregate.sh, gated by the same conditionals,
# so each artifact lands in its proper remote subfolder
# (Ae/, PCA/, Fst/, Froh/, Relatedness/, ROH/bcftools/, het_and_COI/).
#
# Missing files emit a WARNING and the script keeps going — the original
# pipeline would have hard-failed on the same paths.
#
# Usage:
#   OUTPUT_DIR=results_20260421_200003 bash upload_outputs.sh

# Validate before sourcing — common.sh would otherwise default OUTPUT_DIR to a
# new timestamped dir and mkdir it.
if [[ -z "${OUTPUT_DIR:-}" ]]; then
    echo "ERROR: OUTPUT_DIR must be set (e.g. OUTPUT_DIR=results_20260421_200003)" >&2
    exit 1
fi
if [[ ! -d "${OUTPUT_DIR}" ]]; then
    echo "ERROR: OUTPUT_DIR='${OUTPUT_DIR}' is not an existing directory" >&2
    exit 1
fi

# Suppress common.sh's tee-to-run.log so we don't append a fresh upload-only
# session to the original run's log.
export GENDIV_LOG_SETUP=1

source "$(dirname "$0")/genDiversity_common.sh"

# Drop common.sh's ERR trap — every missing file would otherwise emit a
# spurious "Pipeline failed at line ..." line.
trap - ERR

# Cluster module (matches genDiversity_shared.sh). Harmless if `module` is not
# in the environment.
module load rclone 2>/dev/null || true

REMOTE_BASE="remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs"
RCLONE_FLAGS=(--drive-shared-with-me)

# Push a single file. Skip + warn if missing; rclone failures are non-fatal.
push() {
    local src="$1" dest="$2"
    if [[ ! -e "$src" ]]; then
        echo "WARNING: missing (skip): $src" >&2
        return 0
    fi
    rclone -v copy "$src" "${REMOTE_BASE}/${dest}/" "${RCLONE_FLAGS[@]}" \
        || echo "WARNING: upload of $src failed" >&2
}

# Push files in <dir> matching <pattern> (rclone --include glob).
push_glob() {
    local dir="$1" pattern="$2" dest="$3"
    if [[ ! -d "$dir" ]]; then
        echo "WARNING: missing dir (skip): $dir" >&2
        return 0
    fi
    rclone -v copy "$dir" --include "$pattern" "${REMOTE_BASE}/${dest}/" "${RCLONE_FLAGS[@]}" \
        || echo "WARNING: glob upload of $dir / $pattern failed" >&2
}

log "Upload-only run; OUTPUT_DIR=${OUTPUT_DIR}"

# ----------------------------------------------------------------------------
# Per-group uploads (mirrors genDiversity_per_group.sh)
# ----------------------------------------------------------------------------
for rg in wholePop Trotter Pacer; do
    log "Uploading rg=${rg}"
    rg_tag=".${rg}"

    pruned_afreq="${OUTPUT_DIR}/LD_pruned/pruned.${rg}.afreq"
    pca_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune${rg_tag}.pca"
    het_rg_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.het_stats.${rg}"
    freqs_prefix="${OUTPUT_DIR}/divStats/freqs.${rg}"
    roh_RG="${OUTPUT_DIR}/divStats/roh.L3"
    king_prefix="${OUTPUT_DIR}/divStats/filtered.LD_prune.king_gait.${rg}"
    kingkin="${king_prefix}.kin0"
    cross_prefix="${OUTPUT_DIR}/divStats/${rg}.relatedness_correlation"
    froh_prefix="${OUTPUT_DIR}/divStats/coi_Froh_rmDiag_correlation.${rg}"
    froh_cons_prefix="${OUTPUT_DIR}/divStats/coi_Froh_rmDiag_conShare_correlation.${rg}"
    dbl_tag=".${rg}"

    if [[ "$rg" == "wholePop" ]]; then
        subs=("wholePop")
    else
        subs=("$rg" "${rg}_LOW" "${rg}_MEDIUM" "${rg}_HIGH")
    fi

    # §1 per-group afreq
    push "$pruned_afreq" "Ae"

    # §2 PCA + overlays
    push_glob "${OUTPUT_DIR}/divStats" "$(basename "$pca_prefix").eigen*" "PCA"
    push "${OUTPUT_DIR}/divStats/Var_PCs${rg_tag}.jpg" "PCA"
    push "${OUTPUT_DIR}/divStats/pca_plot_BookSize.${rg}.png" "PCA"

    if [[ "$rg" == "wholePop" ]]; then
        # §2 wholePop-only: sex / gait overlays + PC outlier shortlists
        push "${OUTPUT_DIR}/divStats/pca_plot_sex.png" "PCA"
        push "${OUTPUT_DIR}/divStats/pca_plot_Gait.png" "PCA"
        push "${OUTPUT_DIR}/divStats/Trotters_segregating_on_PC2.csv" "PCA"
        push "${OUTPUT_DIR}/divStats/Pacers_cosegregating_withTrotters_on_PC1.csv" "PCA"
        push "${OUTPUT_DIR}/divStats/Trotters_cosegregating_withPacers_on_PC1.csv" "PCA"

        # §1bis Ae (wholePop only)
        push "${OUTPUT_DIR}/divStats/effAllele_stats.txt" "Ae"
        push "${OUTPUT_DIR}/divStats/Figure_Ae_BookSize.tiff" "Ae"

        # §2 whole-pop FST (sex / gait / bookSize concatenated)
        push "${OUTPUT_DIR}/divStats/autosomal.fst.summary" "Fst"
    fi

    # §3 book-size FST within gait (Trotter / Pacer only)
    if [[ "$rg" != "wholePop" ]]; then
        push "${OUTPUT_DIR}/divStats/filtered.LD_prune.fst_bookSize.${rg}.fst.summary" "Fst"
    fi

    # §4 per-group .het (F_SNP)
    push "${het_rg_prefix}.het" "het_and_COI"

    # §4b per-group stratified .het summaries (wholePop + Trotter + Pacer baselines)
    push_glob "${OUTPUT_DIR}/divStats" "filtered.LD_prune.het_stats.${rg}.het.wGait*" "het_and_COI"

    # §5 PCA COI overlay
    push "${OUTPUT_DIR}/divStats/pca_plot_inbreeding.${rg}.png" "PCA"

    # §6 per-group AF table for downstream `bcftools roh --AF-file`
    push "${freqs_prefix}.tab.gz"     "ROH/bcftools"
    push "${freqs_prefix}.tab.gz.tbi" "ROH/bcftools"

    # §6b wholePop-only ROH-vs-.het sanity correlation
    if [[ "$rg" == "wholePop" ]]; then
        push "${OUTPUT_DIR}/divStats/filtered.not_pruned.roh_summary_by_RG_L3.pairplot.png" "ROH/bcftools"
    fi

    # §7 per-base ROH coverage histograms (per sub-group)
    for sub in "${subs[@]}"; do
        push "${roh_RG}.per_base_coverage.${sub}.histo" "ROH/bcftools/freq"
    done

    # §7 consensus ROH summary
    push "${roh_RG}.consensus_${pct}pct.${rg}.summary.txt" "ROH/bcftools"

    # §8 F_ROH summary
    push "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.${rg}.txt" "Froh"

    # §8b wholePop-only F_ROH histogram + high-F_ROH shortlist
    if [[ "$rg" == "wholePop" ]]; then
        push "${OUTPUT_DIR}/divStats/roh_summary_by_RG_L3_Froh.histo" "Froh"
        push "${OUTPUT_DIR}/divStats/roh_high.csv" "Froh"
        push "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait.sumStats.csv" "Froh"
        push "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.sumStats.csv" "Froh"
    fi

    # §9 ROHRM per-cutoff outputs (Inbreeding_Comparison / Pairwise_Differences / matrix png)
    roh_sd_label="${ROH_THRESHOLD_SD%.*}SD"
    for rohrm_mb in $ROH_CUTOFFS; do
        subfolder="${OUTPUT_DIR}/rep_ROHRM/roh_${rohrm_mb%.*}Mb.Threshold_${roh_sd_label}"
        subfolder_tail="$(basename "$subfolder")"
        push "${subfolder}/Robust_Matrix_Comparison_Enhanced.${rg}.png" "Relatedness/${subfolder_tail}"
        push "${subfolder}/Inbreeding_Comparison.${rg}.csv" "Relatedness/${subfolder_tail}"
        push "${subfolder}/Pairwise_Differences.${rg}.csv"  "Relatedness/${subfolder_tail}"
    done

    # §10 per-group KING + IBS
    push "$kingkin"                   "Relatedness"
    push "${kingkin%.kin0}.histo"     "Relatedness"

    # §14 cross-method relatedness correlation
    push "${cross_prefix}.pairplot.png" "Relatedness"

    # §15 froh / froh-cons / doubleAnn correlation plots
    push "${froh_prefix}.pairplot.png"      "Relatedness"
    push "${froh_cons_prefix}.pairplot.png" "Relatedness"
    push "${OUTPUT_DIR}/divStats/correlation_plot_ROH_sh_vs_D_STD_doubleAnn${dbl_tag}.png" "Relatedness"
    push "${OUTPUT_DIR}/divStats/correlation_plot_F_ROH_vs_D_ROH_doubleAnn${dbl_tag}.png" "Relatedness"
    push "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_D_ROH_doubleAnn${dbl_tag}.png" "Relatedness"
    push "${OUTPUT_DIR}/divStats/correlation_plot_F_SNP_vs_F_ROH_doubleAnn${dbl_tag}.png" "Relatedness"
done

# ----------------------------------------------------------------------------
# Cross-group aggregate uploads (mirrors genDiversity_aggregate.sh)
# ----------------------------------------------------------------------------
log "Uploading cross-group aggregate outputs"
push "${OUTPUT_DIR}/divStats/fst_stats.txt" "Fst"

for gp in wholePop twoGait threeBooksize; do
    push "${OUTPUT_DIR}/divStats/Froh_vs_ROHsh_${gp}.png"           "Froh"
    push "${OUTPUT_DIR}/divStats/normalized_ROHsh_${gp}.histogram.png" "Froh"
    push "${OUTPUT_DIR}/divStats/normalized_ROHsh_${gp}.density.png"   "Froh"
    push "${OUTPUT_DIR}/divStats/ROHshared_${gp}.histogram.png"        "Froh"
done

log "Upload-only run finished."
