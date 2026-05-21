#!/usr/bin/env bash
# Per-group SNeP v1.11 wrapper for LD-decay N_e analysis.
#
# SNeP (Barbato et al. 2015; doi:10.3389/fgene.2015.00109) estimates a
# per-generation N_e trajectory from LD decay using the Sved-Feldman
# (1971) approximation. It is the most widely used LD-decay N_e tool in
# 2023-2025 livestock genomics (Manunza 2025 review).
#
# Inputs (CLI):
#   --pruned-prefix <PLINK prefix>     PLINK 1 .bed/.bim/.fam LD-pruned set.
#                                      SNeP bins SNP pairs by recombination
#                                      distance and averages r^2 per bin, so
#                                      LD pruning is recommended (matches
#                                      McGivney 2020 / Manunza 2025 practice).
#   --group         <LABEL:SAMPLES_FILE>   Repeatable; two-column FID/IID.
#   --out-dir       <dir>              <dir>/snep/ is created.
#   --snep-bin      <path>             SNeP binary (default $SNEP_BIN env var
#                                      or tools/SNeP/SNeP1.11 in the repo).
#   --threads       <int>              SNeP -threads (default 8).
#   --maf           <float>            SNeP -maf (default 0.05).
#   --rec-rate-M-bp <float>            Constant recombination rate in M/bp
#                                      (default 1.16e-8 = 1.16 cM/Mb;
#                                      Beeson et al. 2020 EquCab3 mean).
#   --items-threshold <int>            SNeP -itemsTH: minimum SNP pairs per
#                                      distance bin (default 500, SNeP's own
#                                      default). Lower for smoke tests or for
#                                      small populations / sparse panels
#                                      where the default leaves many bins
#                                      empty.
#
# Outputs per group at <dir>/snep/<group>/:
#   * input.<group>.{ped,map}     PLINK PED/MAP fed to SNeP
#   * <group>.NeAll               SNeP per-generation N_e trajectory
#                                 (GenAgo, Ne, d, r^2, ...)
#   * <group>.LDAll               SNeP per-bin LD diagnostic (optional)
#   * <group>.snep.log            SNeP stdout/stderr log
#
# Plus a combined CSV at <dir>/snep/Ne_snep_summary.csv with columns
# Group, Method, Generations_ago, Ne, CI_low_95, CI_high_95 (CIs NA -
# SNeP does not emit per-bin CIs natively).
set -eo pipefail

PRUNED_PREFIX=""
GROUP_SPECS=()
OUT_DIR=""
SNEP_BIN="${SNEP_BIN:-}"
THREADS="8"
MAF="0.05"
REC_RATE="1.16e-8"
ITEMS_TH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pruned-prefix) PRUNED_PREFIX="$2"; shift 2 ;;
        --group)
            GROUP_SPECS+=("$2")
            shift 2
            ;;
        --out-dir)       OUT_DIR="$2"; shift 2 ;;
        --snep-bin)      SNEP_BIN="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --maf)           MAF="$2"; shift 2 ;;
        --rec-rate-M-bp) REC_RATE="$2"; shift 2 ;;
        --items-threshold) ITEMS_TH="$2"; shift 2 ;;
        *) echo "[run_snep] unknown flag: $1" >&2; exit 64 ;;
    esac
done

[[ -z "$PRUNED_PREFIX" ]] && { echo "[run_snep] --pruned-prefix required" >&2; exit 64; }
[[ "${#GROUP_SPECS[@]}" -eq 0 ]] && { echo "[run_snep] at least one --group required" >&2; exit 64; }
[[ -z "$OUT_DIR" ]] && { echo "[run_snep] --out-dir required" >&2; exit 64; }

if [[ -z "$SNEP_BIN" ]]; then
    SNEP_BIN="$(cd "$(dirname "$0")/../.." && pwd)/tools/SNeP/SNeP1.11"
fi
if [[ ! -x "$SNEP_BIN" ]]; then
    echo "[run_snep] ERROR: SNeP binary not found / not executable: $SNEP_BIN" >&2
    echo "[run_snep]   Run scripts/ne/install_snep.sh first." >&2
    exit 2
fi

mkdir -p "$OUT_DIR/snep"
summary_csv="$OUT_DIR/snep/Ne_snep_summary.csv"
echo "Group,Method,Generations_ago,Ne,CI_low_95,CI_high_95" > "$summary_csv"

for spec in "${GROUP_SPECS[@]}"; do
    label="${spec%%:*}"
    samples="${spec#*:}"
    if [[ ! -s "$samples" ]]; then
        echo "[run_snep] WARNING: samples file empty/missing for $label: $samples" >&2
        continue
    fi
    grp_dir="$OUT_DIR/snep/$label"
    mkdir -p "$grp_dir"

    # 1. Subset PLINK and export to PED/MAP. SNeP accepts -ped/-map
    #    universally; -bfile is not exposed by SNeP v1.1x. We use
    #    --output-chr '26' (no chr prefix) because SNeP expects
    #    numeric chromosome codes in the .map file.
    grp_prefix="$grp_dir/input.${label}"
    plink2 --bfile "$PRUNED_PREFIX" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
           --keep "$samples" \
           --export 'ped' \
           --output-chr '26' \
           --out "$grp_prefix" >/dev/null

    # 2. Run SNeP. -recrate is per-bp (1.16e-8 = 1.16 cM/Mb; Beeson
    #    et al. 2020 EquCab3 genome-wide mean). We pick Haldane as the
    #    mapping function (Barbato 2015 convention, also default in the
    #    2023-2025 livestock literature; SNeP requires one of -haldane /
    #    -kosambi / -sved / -svedf). -samplesize 2 applies the diploid
    #    Sved bias correction r2[adj] = r2 - 1/(2n) (Manunza 2025).
    log_file="$grp_dir/snep.${label}.log"
    out_prefix="$grp_dir/${label}"
    snep_args=(
        -ped "$grp_prefix.ped"
        -map "$grp_prefix.map"
        -out "$out_prefix"
        -threads "$THREADS"
        -maf "$MAF"
        -recrate "$REC_RATE"
        -haldane
        -samplesize 2
    )
    [[ -n "$ITEMS_TH" ]] && snep_args+=( -itemsTH "$ITEMS_TH" )
    "$SNEP_BIN" "${snep_args[@]}" > "$log_file" 2>&1 || {
            echo "[run_snep] WARNING: $label SNeP exited non-zero; see $log_file" >&2
            continue
        }

    ne_file="$out_prefix.NeAll"
    if [[ ! -s "$ne_file" ]]; then
        echo "[run_snep] WARNING: $label produced no .NeAll file; see $log_file" >&2
        continue
    fi

    # 3. Standardise the trajectory into the shared CSV. SNeP .NeAll has
    #    a header row (GenAgo Ne d r^2 ...) followed by numeric rows.
    #    We resolve column indices by header name so the parser tolerates
    #    column-order changes between SNeP point releases.
    awk -v g="$label" 'BEGIN{OFS=","}
        NR==1 {
            for (i=1; i<=NF; i++) {
                h = tolower($i); gsub(/[^a-z0-9]/, "", h)
                col[h] = i
            }
            gen_idx = col["genagot"]
            if (!gen_idx) gen_idx = col["genago"]
            if (!gen_idx) gen_idx = col["t"]
            ne_idx  = col["ne"]
            if (!gen_idx || !ne_idx) {
                # Fall back to positional parsing (col1=GenAgo, col2=Ne)
                # if the header is missing or non-standard.
                gen_idx = 1; ne_idx = 2
                if ($1+0 == $1) {
                    # First row is already numeric (no header) — emit it.
                    printf "%s,SNeP,%s,%s,NA,NA\n", g, $gen_idx, $ne_idx
                }
                next
            }
            next
        }
        $1+0 == $1 { printf "%s,SNeP,%s,%s,NA,NA\n", g, $gen_idx, $ne_idx }
    ' "$ne_file" >> "$summary_csv"

    nrows=$(($(wc -l < "$ne_file") - 1))
    echo "[run_snep] $label: SNeP done — $nrows rows in $ne_file"
done

echo "[run_snep] summary written to $summary_csv"
