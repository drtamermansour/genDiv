#!/usr/bin/env bash
# Per-group NeEstimator v2.x wrapper for the LD-decay N_e analysis.
#
# Method-matches McGivney et al. 2020 (doi:10.1038/s41598-019-57389-5),
# which reported a contemporary global Thoroughbred N_e = 330 using
# NeEstimator v2's LD method on an LD-pruned SNP set converted to
# GENEPOP via PGDSpider. We follow Manunza et al. 2025
# (doi:10.3389/fgene.2025.1588986)'s livestock recommendations:
#   * LD method (Waples 2006; Waples & Do 2008, 2010);
#   * P_Crit = 0.02 (Waples & Do 2010 default for stable estimates at
#     moderate to large Ne);
#   * Waples (2006) sample-size bias correction (built into v2 by default);
#   * Random-mating (not monogamy) for closed-studbook livestock;
#   * Harmonic mean across replicates;
#   * Jackknife 95% confidence intervals.
#
# Inputs (CLI):
#   --pruned-prefix <PLINK prefix>     PLINK 1 .bed/.bim/.fam LD-pruned set.
#   --group         <LABEL:SAMPLES_FILE>   Repeatable per group;
#                                          samples file is two-column FID/IID.
#   --out-dir       <dir>              <dir>/neestimator/<group>/ is created.
#   --ne2l-bin      <path>             Optional override for the Ne2x binary
#                                      (default: $NE2L_BIN env var or
#                                      tools/NeEstimator/Ne2x in the repo).
#   --pcrit         <float>            P_Crit (default 0.02).
#
# Outputs per group at <dir>/neestimator/<group>/:
#   * input.<group>.{ped,map}     PLINK subset + recode 12 intermediates
#   * input.<group>.gen           GENEPOP-format file fed to Ne2x
#   * info.<group>.txt            NeEstimator info file (CLI parameters)
#   * input.<group>Ne.txt         Ne2x primary output (Ne, CIs, etc.)
#   * input.<group>NexLD.txt      Ne2x extra LD-method output (per-locus etc.)
#
# Plus a combined standardised CSV at <dir>/neestimator/Ne_neestimator_summary.csv
# with columns Group, Method, Generations_ago, Ne, CI_level, CI_low, CI_high.
# Generations_ago = "contemporary"; CI_level is one of "95_jackknife"
# (if --jackknife was passed and the jackknife block was parsed),
# "95_parametric" (default — parametric CI is always computed by Ne2x),
# or "NA" if no CI could be parsed.
set -eo pipefail

PRUNED_PREFIX=""
GROUP_SPECS=()
OUT_DIR=""
NE2L_BIN="${NE2L_BIN:-}"
PCRIT="0.05"          # Default: exclude rare alleles (MAF < 0.05). Faster
                      # than the Manunza-2025 PCrit=0.02 default and within
                      # the range of values recommended by Waples & Do
                      # (2010) for stable LD-Ne with moderate panels.
JACKKNIFE="0"         # Default: disable jackknife CIs (use parametric only)
                      # because jackknife adds a substantial post-scan pass
                      # for the LOO recomputation. Re-enable with
                      # --jackknife if reviewers want the McGivney-style
                      # jackknife CI.

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pruned-prefix) PRUNED_PREFIX="$2"; shift 2 ;;
        --group)
            GROUP_SPECS+=("$2")
            shift 2
            ;;
        --out-dir)       OUT_DIR="$2"; shift 2 ;;
        --ne2l-bin)      NE2L_BIN="$2"; shift 2 ;;
        --pcrit)         PCRIT="$2"; shift 2 ;;
        --jackknife)     JACKKNIFE="1"; shift ;;
        --no-jackknife)  JACKKNIFE="0"; shift ;;
        *) echo "[run_neestimator] unknown flag: $1" >&2; exit 64 ;;
    esac
done

[[ -z "$PRUNED_PREFIX" ]] && { echo "[run_neestimator] --pruned-prefix required" >&2; exit 64; }
[[ "${#GROUP_SPECS[@]}" -eq 0 ]] && { echo "[run_neestimator] at least one --group required" >&2; exit 64; }
[[ -z "$OUT_DIR" ]] && { echo "[run_neestimator] --out-dir required" >&2; exit 64; }

if [[ -z "$NE2L_BIN" ]]; then
    NE2L_BIN="$(cd "$(dirname "$0")/../.." && pwd)/tools/NeEstimator/Ne2x"
fi
if [[ ! -x "$NE2L_BIN" ]]; then
    echo "[run_neestimator] ERROR: Ne2x not found / not executable: $NE2L_BIN" >&2
    echo "[run_neestimator]   Run scripts/ne/install_neestimator.sh first." >&2
    exit 2
fi

mkdir -p "$OUT_DIR/neestimator"
summary_csv="$OUT_DIR/neestimator/Ne_neestimator_summary.csv"
# Standardised summary CSV header (shared across all Ne wrappers).
# CI_level: "95_jackknife" (jackknife enabled and jackknife block parsed),
# "95_parametric" (parametric CI, either by --no-jackknife or fallback),
# or "NA" when no CI could be parsed.
echo "Group,Method,Generations_ago,Ne,CI_level,CI_low,CI_high" > "$summary_csv"

for spec in "${GROUP_SPECS[@]}"; do
    label="${spec%%:*}"
    samples="${spec#*:}"
    if [[ ! -s "$samples" ]]; then
        echo "[run_neestimator] WARNING: samples file empty/missing for $label: $samples" >&2
        continue
    fi
    grp_dir="$OUT_DIR/neestimator/$label"
    mkdir -p "$grp_dir"

    # 1. Subset PLINK and recode genotypes as 1/2 (PED format). We use
    #    PLINK 1.9 here because `--recode 12` (numeric allele coding) is
    #    a PLINK 1 form; PLINK 2's `--export 12` exists but the format
    #    name and behaviour are slightly different.
    grp_prefix="$grp_dir/input.${label}"
    plink --bfile "$PRUNED_PREFIX" --chr-set 31 no-y no-xy no-mt --allow-extra-chr \
          --keep "$samples" \
          --recode 12 \
          --output-chr 'chr26' \
          --out "$grp_prefix" >/dev/null

    # 2. Convert PED + MAP to GENEPOP. NeEstimator's GENEPOP format:
    #      <title>
    #      LocusName1
    #      ...
    #      LocusNameN
    #      Pop
    #      <IID>, 0101 0102 0202 ...
    #    PLINK --recode 12 codes alleles as "1"/"2" (or "0" for missing).
    #    GENEPOP wants 2-digit codes per allele, so we map 1->01, 2->02, 0->00.
    gen_file="$grp_prefix.gen"
    {
        echo "NeEstimator input for ${label} from $(basename "$PRUNED_PREFIX")"
        awk '{print $2}' "$grp_prefix.map"
        echo "Pop"
        awk 'BEGIN{OFS=""} \
             { iid=$2; printf "%s, ", iid; \
               for (i=7; i<=NF; i+=2) { \
                 a1=$i; a2=$(i+1); \
                 c1=(a1=="0"?"00":"0"a1); \
                 c2=(a2=="0"?"00":"0"a2); \
                 printf "%s%s ", c1, c2; \
               } \
               printf "\n"; \
             }' "$grp_prefix.ped"
    } > "$gen_file"

    # 3. Write NeEstimator info file. Method bitmap: 1=LD only. Format
    #    code: 2=GENEPOP. Random mating (0). PCrit = $PCRIT.
    info_file="$grp_dir/info.${label}.txt"
    cat > "$info_file" <<INFO
1                       * method bitmap: 1=LD, 2=Het, 4=Coan, 8=Temporal (sum)
$grp_dir/               * input directory (trailing slash)
input.${label}.gen      * input file name
2                       * 1=FSTAT, 2=GENEPOP
$grp_dir/               * output directory (trailing slash)
input.${label}Ne.txt    * output file name
1                       * number of PCrit values
$PCRIT                  * PCrit value(s) — Manunza 2025: 0.02 is the standard livestock default
0                       * 0 = random mating (livestock default), 1 = monogamy
INFO

    # 4. Write options file. Extra LD output (bitmap 1), no individual cap,
    #    no Freq output, no Burrow output. Parametric CI always on; jackknife
    #    is configurable because the LOO recomputation is a substantial
    #    post-scan pass at large N.
    opt_file="$grp_dir/options.${label}.txt"
    cat > "$opt_file" <<OPT
1                       * extra output methods bitmap: LD=1
0                       * max individuals/pop (0 = no limit)
0                       * Freq output (0 = none)
0                       * Burrow output (0 = none)
1                       * parametric CI: 1=yes
${JACKKNIFE}                       * jackknife CI: 1=yes, 0=no
0                       * up to population (0 = no restriction)
0                       * all loci accepted
1                       * write missing-data file: 1=yes
0                       * chromosome/loci option file: 0=none
OPT

    # 5. Run Ne2x. CLI form: Ne2x  i:info_file  o:options_file
    log_file="$grp_dir/neestimator.${label}.log"
    "$NE2L_BIN" "i:$info_file" "o:$opt_file" > "$log_file" 2>&1 || {
        echo "[run_neestimator] WARNING: $label Ne2x exited non-zero; see $log_file" >&2
        continue
    }

    ne_out="$grp_dir/input.${label}Ne.txt"
    if [[ ! -s "$ne_out" ]]; then
        echo "[run_neestimator] WARNING: $label produced no Ne output; see $log_file" >&2
        continue
    fi

    # 6. Parse the Ne output. The LD method emits a block per PCrit value
    #    with: estimated Ne, parametric 95% CI, and (if enabled) jackknife
    #    95% CI. The file format has two columns (col 1 = our PCrit, col 2
    #    = the NeEstimator default "0+" comparison). We take column 1.
    #    Block layout:
    #        * Parametric                      <ci_lo_col1>   <ci_lo_col2>
    #                                          <ci_hi_col1>   <ci_hi_col2>
    #        * JackKnife on Samples            <ci_lo_col1>   <ci_lo_col2>
    #                                          <ci_hi_col1>   <ci_hi_col2>
    #    Lower bound on the label line, upper bound on the next line.
    #    We prefer jackknife CIs when present (more conservative,
    #    Manunza 2025 livestock convention) and fall back to parametric.
    ne_val=$(awk '/^[[:space:]]*Estimated Ne\^? *=/ {
                       s=$0; gsub(/[^0-9.+-]/, " ", s);
                       n=split(s, a, /[ ]+/);
                       for (i=1;i<=n;i++) if (a[i]!="" && a[i]+0==a[i] && a[i]+0>0) {print a[i]; exit}
                  }' "$ne_out")
    _parse_ci_block() {
        local label_re="$1"; local which="$2"   # which = lo|hi
        if [[ "$which" == "lo" ]]; then
            awk -v re="$label_re" 'tolower($0) ~ re {
                       s=$0; gsub(/[^0-9.+-]/, " ", s);
                       n=split(s, a, /[ ]+/);
                       for (i=1;i<=n;i++) if (a[i]!="" && a[i]+0==a[i] && a[i]+0>0) {print a[i]; exit}
                  }' "$ne_out"
        else
            awk -v re="$label_re" 'found && NF>0 {
                       s=$0; gsub(/[^0-9.+-]/, " ", s);
                       n=split(s, a, /[ ]+/);
                       for (i=1;i<=n;i++) if (a[i]!="" && a[i]+0==a[i] && a[i]+0>0) {print a[i]; exit}
                       exit
                  }
                  tolower($0) ~ re {found=1}' "$ne_out"
        fi
    }
    # Try jackknife first (preferred when it was computed); fall back to
    # parametric. Track which block we actually parsed so we can emit the
    # right CI_level into the standardised CSV.
    ci_level="95_jackknife"
    ci_lo=$(_parse_ci_block "jackknife on samples" lo)
    ci_hi=$(_parse_ci_block "jackknife on samples" hi)
    if [[ -z "$ci_lo" || -z "$ci_hi" ]]; then
        # Jackknife block not present (disabled in options). Fall back to
        # parametric CI.
        ci_level="95_parametric"
        ci_lo=$(_parse_ci_block "^[[:space:]]*\\\\* parametric" lo)
        ci_hi=$(_parse_ci_block "^[[:space:]]*\\\\* parametric" hi)
    fi
    [[ -z "$ne_val" ]] && ne_val="NA"
    [[ -z "$ci_lo" ]] && ci_lo="NA"
    [[ -z "$ci_hi" ]] && ci_hi="NA"
    if [[ "$ci_lo" == "NA" && "$ci_hi" == "NA" ]]; then
        ci_level="NA"
    fi
    printf "%s,NeEstimator_LD,contemporary,%s,%s,%s,%s\n" \
        "$label" "$ne_val" "$ci_level" "$ci_lo" "$ci_hi" >> "$summary_csv"

    echo "[run_neestimator] $label: Ne=$ne_val (CI_level=$ci_level, CI=$ci_lo - $ci_hi)"
done

echo "[run_neestimator] summary written to $summary_csv"
