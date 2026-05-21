#!/usr/bin/env bash
# Per-group currentNe2 wrapper for the contemporary single-point Ne
# analysis. currentNe2 (Santiago, Köpke & Caballero 2025;
# doi:10.1038/s41467-025-61378-w) estimates contemporary Ne from LD
# between (mostly unlinked) loci. It is the single-point contemporary
# analogue to GONE2 and is recommended by the Santiago 2025 paper as a
# faster alternative to NeEstimator v2 for the same use case.
#
# Inputs (CLI):
#   --pruned-prefix <PLINK prefix>     PLINK 1 .bed/.bim/.fam LD-pruned set.
#                                      currentNe2 benefits from a panel of
#                                      mostly-unlinked SNPs, so LD pruning
#                                      is appropriate here (same input we
#                                      feed to NeEstimator).
#   --group         <LABEL:SAMPLES_FILE>   Repeatable; two-column FID/IID.
#   --out-dir       <dir>              <dir>/currentne2/ is created.
#   --currentne2-bin <path>            Optional; default $CURRENTNE2_BIN env
#                                      or tools/currentNe2/currentne2.
#   --rec-rate-cm-mb <float>           cM/Mb conversion (default 1.16;
#                                      Beeson et al. 2020 EquCab3 mean).
#                                      Used by currentNe2 -r to convert
#                                      physical to genetic positions.
#   --threads       <int>              currentNe2 -t (default 8).
#   --extra-args    <string>           Pass-through string of additional
#                                      currentNe2 flags (e.g. "-k -1"
#                                      for sample-average sib correction).
#
# Outputs per group at <dir>/currentne2/<group>/:
#   * input.<group>.{ped,map}             PLINK subset fed to currentNe2
#   * input.<group>_currentNe2_OUTPUT.txt currentNe2 primary output
#   * currentne2.<group>.log              stdout/stderr log
#
# Plus a combined CSV at <dir>/currentne2/Ne_currentne2_summary.csv with
# columns Group, Method, Generations_ago, Ne, CI_low_95, CI_high_95.
set -eo pipefail

PRUNED_PREFIX=""
GROUP_SPECS=()
OUT_DIR=""
CURRENTNE2_BIN="${CURRENTNE2_BIN:-}"
REC_RATE="1.16"
THREADS="8"
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pruned-prefix)   PRUNED_PREFIX="$2"; shift 2 ;;
        --group)
            GROUP_SPECS+=("$2"); shift 2
            ;;
        --out-dir)         OUT_DIR="$2"; shift 2 ;;
        --currentne2-bin)  CURRENTNE2_BIN="$2"; shift 2 ;;
        --rec-rate-cm-mb)  REC_RATE="$2"; shift 2 ;;
        --threads)         THREADS="$2"; shift 2 ;;
        --extra-args)      EXTRA_ARGS="$2"; shift 2 ;;
        *) echo "[run_currentne2] unknown flag: $1" >&2; exit 64 ;;
    esac
done

[[ -z "$PRUNED_PREFIX" ]] && { echo "[run_currentne2] --pruned-prefix required" >&2; exit 64; }
[[ "${#GROUP_SPECS[@]}" -eq 0 ]] && { echo "[run_currentne2] at least one --group required" >&2; exit 64; }
[[ -z "$OUT_DIR" ]] && { echo "[run_currentne2] --out-dir required" >&2; exit 64; }

if [[ -z "$CURRENTNE2_BIN" ]]; then
    CURRENTNE2_BIN="$(cd "$(dirname "$0")/../.." && pwd)/tools/currentNe2/currentne2"
fi
if [[ ! -x "$CURRENTNE2_BIN" ]]; then
    echo "[run_currentne2] ERROR: currentNe2 binary not found / not executable: $CURRENTNE2_BIN" >&2
    echo "[run_currentne2]   Run scripts/ne/install_currentne2.sh first." >&2
    exit 2
fi

mkdir -p "$OUT_DIR/currentne2"
summary_csv="$OUT_DIR/currentne2/Ne_currentne2_summary.csv"
echo "Group,Method,Generations_ago,Ne,CI_low_95,CI_high_95" > "$summary_csv"

for spec in "${GROUP_SPECS[@]}"; do
    label="${spec%%:*}"
    samples="${spec#*:}"
    if [[ ! -s "$samples" ]]; then
        echo "[run_currentne2] WARNING: samples file empty/missing for $label: $samples" >&2
        continue
    fi
    grp_dir="$OUT_DIR/currentne2/$label"
    mkdir -p "$grp_dir"

    # 1. Subset PLINK to this group's samples and export to PED/MAP.
    grp_prefix="$grp_dir/input.${label}"
    plink2 --bfile "$PRUNED_PREFIX" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
           --keep "$samples" \
           --export 'ped' \
           --output-chr 'chr26' \
           --out "$grp_prefix" >/dev/null

    # 2. Run currentNe2. The PED form takes a constant cM/Mb rec rate via
    #    -r; we use Beeson 2020's EquCab3 genome-wide mean of 1.16 cM/Mb.
    #    -t selects threads (OpenMP-parallel). Output filename uses the
    #    suffix _currentNe2_OUTPUT.txt automatically.
    log_file="$grp_dir/currentne2.${label}.log"
    cmd=( "$CURRENTNE2_BIN" -r "$REC_RATE" -t "$THREADS" )
    if [[ -n "$EXTRA_ARGS" ]]; then
        # shellcheck disable=SC2206
        cmd+=( $EXTRA_ARGS )
    fi
    cmd+=( "$grp_prefix.ped" )
    echo "[run_currentne2] $label: ${cmd[*]}" > "$log_file"
    "${cmd[@]}" >> "$log_file" 2>&1 || {
        echo "[run_currentne2] WARNING: $label currentNe2 exited non-zero; see $log_file" >&2
        continue
    }

    ne_out="${grp_prefix}_currentNe2_OUTPUT.txt"
    if [[ ! -s "$ne_out" ]]; then
        echo "[run_currentne2] WARNING: $label produced no _currentNe2_OUTPUT.txt; see $log_file" >&2
        continue
    fi

    # 3. Parse the currentNe2 output. The file emits two estimates: a
    #    whole-genome integration over all SNP pairs, and an estimate
    #    using only inter-chromosomal pairs (c = 0.5, truly unlinked).
    #    Following Santiago 2025, we report the inter-chromosomal value
    #    as the primary contemporary Ne (purer because no within-chromosome
    #    LD-decay contamination). CIs reported by currentNe2 are 50% and
    #    90% — we take the 90% values; they live in the CI_low_95 /
    #    CI_high_95 columns of the standardised summary CSV (column names
    #    are kept for cross-tool consistency; see CI_level column).
    #    The file layout has each label on a "#" line followed by the
    #    numeric value on the next non-comment line.
    parse_block_value() {
        # $1 = section start regex, $2 = key regex inside section
        awk -v sect="$1" -v key="$2" '
            $0 ~ sect { in_sect=1; next }
            in_sect && $0 ~ key { want=1; next }
            in_sect && want && /^[[:space:]]*[0-9.+-]+[[:space:]]*$/ {
                gsub(/[[:space:]]/, "", $0); print $0; exit
            }
        ' "$ne_out"
    }
    ne_val=$(parse_block_value "Ne estimation based only on LD between chromosomes" "Ne point estimate")
    ci_lo=$(parse_block_value  "Ne estimation based only on LD between chromosomes" "Lower limit of 90% CI")
    ci_hi=$(parse_block_value  "Ne estimation based only on LD between chromosomes" "Upper limit of 90% CI")
    # Also capture the whole-genome estimate for diagnostic reporting.
    ne_val_whole=$(parse_block_value "Ne estimation by integration over the whole genome" "Ne point estimate")

    [[ -z "$ne_val" ]] && ne_val="NA"
    [[ -z "$ci_lo" ]]  && ci_lo="NA"
    [[ -z "$ci_hi" ]]  && ci_hi="NA"
    [[ -z "$ne_val_whole" ]] && ne_val_whole="NA"
    printf "%s,currentNe2,contemporary,%s,%s,%s\n" \
        "$label" "$ne_val" "$ci_lo" "$ci_hi" >> "$summary_csv"

    echo "[run_currentne2] $label: Ne_interchrom=$ne_val (90% CI: $ci_lo - $ci_hi); Ne_wholegenome=$ne_val_whole"
done

echo "[run_currentne2] summary written to $summary_csv"
