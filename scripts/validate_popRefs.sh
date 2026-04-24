#!/usr/bin/env bash
# Validator for the per-group population-reference files.
#
# Runs in two modes:
#   --mode upstream    Check the pipeline's OUTPUT_DIR, where files sit in
#                      divStats/, LD_pruned/, rep_ROHRM/roh_1Mb.Threshold_3SD/,
#                      and preprocess/.
#   --mode popFiles    Check a flat popFiles/ directory as produced by GPA's
#                      create_popFiles.sh — every file has the same basename
#                      but lives directly under --root.
#
# In either mode it iterates over rg ∈ {wholePop, Trotter, Pacer}, asserts each
# expected file exists, is non-empty, has the expected header, and (where
# applicable) has a row count within sanity bounds. Exits 0 if every check
# passes; exits 1 with a terse failure report otherwise.
#
# GPA-side usage:
#   1. Copy this file to the GPA repo (anywhere on PATH or in scripts/).
#   2. Run it at the end of create_popFiles.sh:
#        bash scripts/validate_popRefs.sh --mode popFiles --root "$popFiles_dir"
#   3. If you need to add or drop files checked, edit the FILE_SPECS table
#      below — everything else is parameterised.

set -eo pipefail

mode=""
root=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)    mode="$2"; shift 2 ;;
        --root)    root="$2"; shift 2 ;;
        -h|--help) sed -n '2,23p' "$0"; exit 0 ;;
        *)         echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done
[[ -z "$mode" || -z "$root" ]] && { echo "usage: $0 --mode upstream|popFiles --root DIR" >&2; exit 2; }
[[ ! -d "$root" ]] && { echo "ERROR: --root is not a directory: $root" >&2; exit 2; }
case "$mode" in upstream|popFiles) : ;; *) echo "ERROR: --mode must be upstream or popFiles" >&2; exit 2 ;; esac

##############################################################################
# FILE_SPECS — single source of truth for what we validate.
# Each entry: logical_key | upstream_subdir | filename_template | header_regex | min_rows_formula
#   filename_template uses %RG% as the group placeholder.
#   header_regex is extended-regex (grep -E) that the first non-comment line
#     must match. Use '' for files with no header (BED files).
#   min_rows_formula is an arithmetic expression over NGRP (group size) that
#     the file must meet or exceed; use 'any' to skip the row-count check.
#     Uses awk arithmetic; e.g. "NGRP+1" for one-per-sample files, or
#     "(NGRP*(NGRP-1))/2+1" for pairwise files.
# The upstream_subdir column is ignored in --mode popFiles (everything is flat).
##############################################################################
FILE_SPECS=(
    "afreq|LD_pruned|pruned.%RG%.afreq|^#CHROM[[:space:]]+ID[[:space:]]+REF[[:space:]]+ALT|any"
    "freqs|divStats|freqs.%RG%.tab.gz||any"
    "freqs_tbi|divStats|freqs.%RG%.tab.gz.tbi||any"
    "het|divStats|filtered.LD_prune.het_stats.%RG%.het|^#?FID[[:space:]]+IID[[:space:]]+O\\(HOM\\)|NGRP+1"
    "froh|divStats|roh_summary_by_RG_L3_Froh.%RG%.txt|^IID[[:space:]]+NSEG[[:space:]]+KB[[:space:]]+KBAVG[[:space:]]+F_ROH$|NGRP+1"
    "inbreeding|rep_ROHRM/roh_1Mb.Threshold_3SD|Inbreeding_Comparison.%RG%.csv|^IID,D_STD,D_ROH,Phenotype$|NGRP+1"
    "pairwise|rep_ROHRM/roh_1Mb.Threshold_3SD|Pairwise_Differences.%RG%.csv|^ID1,ID2,Pheno1,Pheno2,Kinship_Std,Kinship_ROH,Difference|(NGRP*(NGRP-1))/2+1"
    "consensus|divStats|roh.L3.consensus_25pct.merged.%RG%.smoothed.bed||any"
    "rohsh|divStats|roh.L3.perSample_intersect_%RG%_consensus_25pct.summary.txt|^IID[[:space:]]+Total_ROH_in_Consensus_region\\(bp\\)[[:space:]]+Percent_of_Consensus_ROH$|NGRP+1"
)

# sample_groups.tsv is a global (single-file) reference. Checked once, not per rg.
GLOBAL_SPECS=(
    "sample_groups|preprocess|sample_groups.tsv|^IID[[:space:]]+group$|any"
)

##############################################################################

pass=0; fail=0
fail_list=()

resolve_path() {
    local subdir="$1" filename="$2"
    if [[ "$mode" == "upstream" ]]; then
        printf '%s' "${root}/${subdir}/${filename}"
    else
        printf '%s' "${root}/${filename}"
    fi
}

group_size() {
    # Count samples in preprocess/samples.${rg}.txt relative to --root. In
    # popFiles mode the samples file may or may not be present; if missing we
    # fall back to counting rows of whatever per-sample file we can find.
    local rg="$1"
    local upstream_samples="${root}/preprocess/samples.${rg}.txt"
    local popfiles_samples="${root}/samples.${rg}.txt"
    if [[ -f "$upstream_samples" ]]; then
        wc -l < "$upstream_samples" | tr -d ' '
    elif [[ -f "$popfiles_samples" ]]; then
        wc -l < "$popfiles_samples" | tr -d ' '
    else
        echo 0
    fi
}

eval_rows_formula() {
    local formula="$1" ngrp="$2"
    [[ "$formula" == "any" ]] && { echo -1; return; }
    awk -v N="$ngrp" "BEGIN{print ${formula//NGRP/N}}"
}

check_file() {
    local label="$1" path="$2" header_rx="$3" min_rows="$4"
    if [[ ! -f "$path" ]]; then
        fail_list+=("MISSING  $label  ($path)")
        fail=$((fail+1))
        return
    fi
    if [[ ! -s "$path" ]]; then
        fail_list+=("EMPTY    $label  ($path)")
        fail=$((fail+1))
        return
    fi
    if [[ -n "$header_rx" ]]; then
        # Scan the first 10 lines for any match. Handles both PLINK2 files
        # where "#CHROM" is itself the header and files that emit a few
        # "#"-commented lines before the real header (e.g. sample_groups.tsv).
        if ! head -n 10 "$path" | grep -Eq "$header_rx"; then
            local first
            first=$(head -n1 "$path")
            fail_list+=("BAD_HEADER  $label  expected /$header_rx/  (first line: '$first')")
            fail=$((fail+1))
            return
        fi
    fi
    if [[ "$min_rows" != "-1" ]]; then
        local actual
        actual=$(wc -l < "$path" | tr -d ' ')
        if [[ "$actual" -lt "$min_rows" ]]; then
            fail_list+=("TOO_FEW_ROWS  $label  got=$actual expected>=$min_rows  ($path)")
            fail=$((fail+1))
            return
        fi
    fi
    pass=$((pass+1))
}

echo "Validating popRefs (--mode $mode --root $root)"
echo

# Per-group files
for rg in wholePop Trotter Pacer; do
    ngrp=$(group_size "$rg")
    echo "-- rg=$rg  (n=$ngrp)"
    for spec in "${FILE_SPECS[@]}"; do
        IFS='|' read -r key subdir tmpl hdr_rx rows_f <<<"$spec"
        filename="${tmpl//%RG%/$rg}"
        path=$(resolve_path "$subdir" "$filename")
        min_rows=$(eval_rows_formula "$rows_f" "$ngrp")
        check_file "${key}/${rg}" "$path" "$hdr_rx" "$min_rows"
    done
done

# Global files
echo "-- global"
for spec in "${GLOBAL_SPECS[@]}"; do
    IFS='|' read -r key subdir tmpl hdr_rx rows_f <<<"$spec"
    path=$(resolve_path "$subdir" "$tmpl")
    min_rows=$(eval_rows_formula "$rows_f" "0")
    check_file "$key" "$path" "$hdr_rx" "$min_rows"
done

echo
echo "==================================================================="
echo "Passed: $pass"
echo "Failed: $fail"
if [[ "$fail" -gt 0 ]]; then
    printf '  %s\n' "${fail_list[@]}"
    exit 1
fi
echo "All popRefs present and well-shaped."
