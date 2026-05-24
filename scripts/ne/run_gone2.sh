#!/usr/bin/env bash
# Per-group GONE2 wrapper for the LD-decay N_e analysis.
#
# GONE2 (Santiago, Köpke & Caballero 2025; doi:10.1038/s41467-025-61378-w)
# estimates a per-generation N_e trajectory from a single SNP sample,
# using the genome-wide LD spectrum across recombination-distance bins.
#
# Inputs (CLI):
#   --unpruned-prefix <PLINK prefix>   PLINK 1 .bed/.bim/.fam set after QC
#                                      but before LD pruning. GONE2 uses the
#                                      full LD spectrum so pruning would
#                                      discard the signal it estimates from.
#   --group           <LABEL:SAMPLES_FILE>
#                                      Repeatable per group; each samples
#                                      file is the standard PLINK FID/IID
#                                      two-column format.
#   --out-dir         <dir>            Output directory; <dir>/gone2/ is
#                                      created.
#   --gone2-bin       <path>           Path to the GONE2 binary (default:
#                                      env var $GONE2_BIN or
#                                      tools/GONE2/gone2 in the repo).
#   --rec-rate-cm-mb  <float>          cM/Mb conversion (default 1.16;
#                                      Beeson et al. 2020 EquCab3 mean).
#   --threads         <int>            GONE2 -t flag (default 8).
#   --maf             <float>          GONE2 -M flag (default 0.05).
#
# Outputs:
#   <out-dir>/gone2/<group>/           Raw GONE2 outputs per group
#                                        (*_GONE_Ne, *_GONE_d2, *_GONE_STATS).
#   <out-dir>/gone2/Ne_gone2_summary.csv
#                                      Standardised per-group N_e trajectory
#                                      with columns Group, Method,
#                                      Generations_ago, Ne, CI_level,
#                                      CI_low, CI_high. (No native CIs
#                                      from GONE2; CI_level/CI_low/CI_high
#                                      are all emitted as NA.)
set -eo pipefail

UNPRUNED_PREFIX=""
GROUP_SPECS=()
OUT_DIR=""
REC_RATE="1.16"
THREADS="8"
MAF="0.05"
GONE2_BIN="${GONE2_BIN:-}"
METAPOPULATION=0          # GONE2 -x: assume sample is from a metapopulation
                          # with subpopulations of equal size. Default off
                          # (panmixia); pass --metapopulation to enable.
OUTPUT_SUBDIR="gone2"     # Tool sub-directory under --out-dir. Override to
                          # e.g. gone2_x when running a sensitivity analysis
                          # with -x without clobbering the panmixia outputs.
METHOD_LABEL="GONE2"      # Label used in the standardised summary CSV's
                          # Method column. Override to e.g. GONE2_x so the
                          # downstream plotter / summary can distinguish.

while [[ $# -gt 0 ]]; do
    case "$1" in
        --unpruned-prefix) UNPRUNED_PREFIX="$2"; shift 2 ;;
        --group)
            GROUP_SPECS+=("$2")
            shift 2
            ;;
        --out-dir)         OUT_DIR="$2"; shift 2 ;;
        --gone2-bin)       GONE2_BIN="$2"; shift 2 ;;
        --rec-rate-cm-mb)  REC_RATE="$2"; shift 2 ;;
        --threads)         THREADS="$2"; shift 2 ;;
        --maf)             MAF="$2"; shift 2 ;;
        --metapopulation)  METAPOPULATION=1; shift ;;
        --output-subdir)   OUTPUT_SUBDIR="$2"; shift 2 ;;
        --method-label)    METHOD_LABEL="$2"; shift 2 ;;
        *) echo "[run_gone2] unknown flag: $1" >&2; exit 64 ;;
    esac
done

[[ -z "$UNPRUNED_PREFIX" ]] && { echo "[run_gone2] --unpruned-prefix required" >&2; exit 64; }
[[ "${#GROUP_SPECS[@]}" -eq 0 ]] && { echo "[run_gone2] at least one --group LABEL:SAMPLES required" >&2; exit 64; }
[[ -z "$OUT_DIR" ]] && { echo "[run_gone2] --out-dir required" >&2; exit 64; }

if [[ -z "$GONE2_BIN" ]]; then
    GONE2_BIN="$(cd "$(dirname "$0")/../.." && pwd)/tools/GONE2/gone2"
fi
if [[ ! -x "$GONE2_BIN" ]]; then
    echo "[run_gone2] ERROR: GONE2 binary not found / not executable: $GONE2_BIN" >&2
    echo "[run_gone2]   Run scripts/ne/install_gone2.sh first." >&2
    exit 2
fi

mkdir -p "$OUT_DIR/$OUTPUT_SUBDIR"
summary_csv="$OUT_DIR/$OUTPUT_SUBDIR/Ne_${OUTPUT_SUBDIR}_summary.csv"
# Standardised summary CSV header (shared across all Ne wrappers in this
# pipeline). CI_level records what level the CI bounds correspond to and
# how they were computed, so a downstream reader can distinguish
# parametric-95% (NeEstimator), jackknife-95% (NeEstimator with
# --jackknife), 90% (currentNe2), or NA (GONE2 / GONE2_x / SNeP, which do
# not emit CIs natively).
echo "Group,Method,Generations_ago,Ne,CI_level,CI_low,CI_high" > "$summary_csv"

for spec in "${GROUP_SPECS[@]}"; do
    label="${spec%%:*}"
    samples="${spec#*:}"
    if [[ ! -s "$samples" ]]; then
        echo "[run_gone2] WARNING: samples file empty/missing for $label: $samples" >&2
        continue
    fi
    grp_dir="$OUT_DIR/$OUTPUT_SUBDIR/$label"
    mkdir -p "$grp_dir"

    # 1. Subset PLINK to this group's samples; convert to PED/MAP for GONE2.
    grp_prefix="$grp_dir/input.${label}"
    plink2 --bfile "$UNPRUNED_PREFIX" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
           --keep "$samples" \
           --recode 'ped' \
           --output-chr 'chr26' \
           --out "$grp_prefix" >/dev/null

    # 2. Run GONE2. The .map sidecar must sit next to the .ped; PLINK already
    #    placed it there. -r sets a constant cM/Mb rate (we use the EquCab3
    #    genome-wide mean of 1.16 cM/Mb per Beeson et al. 2020); -g 0 forces
    #    unphased-diploid mode (default); -M screens MAF; -t parallel threads.
    # -x: assumes the sample is from a metapopulation with subpopulations of
    #    equal size (Santiago 2025). Off by default (panmixia).
    log_file="$grp_dir/gone2.${label}.log"
    gone2_args=( -g 0 -r "$REC_RATE" -t "$THREADS" -M "$MAF" -o "$grp_dir/${label}" )
    [[ "$METAPOPULATION" == "1" ]] && gone2_args+=( -x )
    gone2_args+=( "$grp_prefix.ped" )
    "$GONE2_BIN" "${gone2_args[@]}" > "$log_file" 2>&1

    # GONE2 writes output as either ${label}_GONE2_Ne (panmixia / default)
    # or ${label}_GONE2_Ne_mix (-x metapopulation mode). The two formats
    # differ — see standardisation block below.
    if [[ "$METAPOPULATION" == "1" ]]; then
        ne_file="$grp_dir/${label}_GONE2_Ne_mix"
    else
        ne_file="$grp_dir/${label}_GONE2_Ne"
    fi
    if [[ ! -s "$ne_file" ]]; then
        echo "[run_gone2] WARNING: $label produced no Ne file ($ne_file); check log:" >&2
        echo "[run_gone2]   $log_file" >&2
        continue
    fi

    # Standardise the trajectory into the shared CSV.
    #   *_GONE2_Ne     (panmixia):       2 cols — Generation, Ne_diploids.
    #   *_GONE2_Ne_mix (metapopulation): 5 cols — Rec_rate_bin, generation,
    #     N_T_metapop, Ne_metapop, d2_metapop, preceded by a multi-line
    #     metadata block (Number of subpopulations / FST / migration rate).
    # We pick the right columns by format.
    if [[ "$METAPOPULATION" == "1" ]]; then
        awk -v g="$label" -v m="$METHOD_LABEL" 'BEGIN{OFS=","} \
            /^[0-9.]/ && NF >= 4 { printf "%s,%s,%s,%s,NA,NA,NA\n", g, m, $2, $4 }' \
            "$ne_file" >> "$summary_csv"
    else
        awk -v g="$label" -v m="$METHOD_LABEL" 'BEGIN{OFS=","} \
            /^[0-9]/ { printf "%s,%s,%s,%s,NA,NA,NA\n", g, m, $1, $2 }' \
            "$ne_file" >> "$summary_csv"
    fi

    echo "[run_gone2] $label: GONE2 done — $(wc -l < "$ne_file" | tr -d " ") rows in $ne_file"
done

echo "[run_gone2] summary written to $summary_csv"
