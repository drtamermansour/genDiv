# Shared configuration, helper functions, and logging setup for the genDiv pipeline.
# Sourced (not executed) by genDiversity.sh and its subscripts.

# ============================================================
# CONFIG — edit these values to adapt the pipeline
# ============================================================
# Paths (relative to the parent directory of this repo)
equCab3_map="$(pwd)/../Equine80select_remapper/results_E80selv2_to_equCab3noAlt_genDiv/qc/Equine80select_v2_1_HTS_20143333_B1_UCD_allele_map_equCab3noAlt.tsv"
ref="../Horse_parentage_SNPs/equCab3/download/equCab3.fa"
reference_fai="../Horse_parentage_SNPs/equCab3/equCab3_genome.fa.fai"
GDRIVE_BASE="remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs"

# Threads and reproducibility
nthreads=10
BEAGLE_SEED=12345       # Fixed seed for reproducible phasing

# Quality control thresholds
KING_CUTOFF=0.20        # Relatedness cutoff for removing first-degree relatives
HWE_PVAL=1e-6           # Hardy-Weinberg equilibrium p-value threshold
GENO_MISS=0.05          # Variant and sample missingness threshold
MAF=0.01                # Minor allele frequency threshold
SEX_MAX_FEMALE_XF=0.2   # Max X-chr inbreeding coefficient to call female
SEX_MIN_MALE_XF=0.8     # Min X-chr inbreeding coefficient to call male
PAR_END_BP=2063653       # End position of pseudoautosomal region on X chr (EquCab3)

# LD pruning
LD_WINDOW_KB=100        # Sliding window size (kb)
LD_R2=0.8               # R² threshold

# ROH and consensus parameters
ROH_THRESHOLD_SD=3.0    # Window filter: Mean - (k*SD) SNP count cutoff
pct=25                  # Minimum % of samples in ROH to define consensus
CONSENSUS_MIN_MB=0.5    # Minimum consensus ROH region size (Mb)
PRIMARY_ROH_MB=1.0      # Primary ROH window cutoff used in downstream analyses
ROH_CUTOFFS="1.0 5.0 10.0"   # All ROH window cutoffs to evaluate

# Directory layout (relative to the working directory)
OUTPUT_DIR="${OUTPUT_DIR:-results_$(date +%Y%m%d_%H%M%S)}" # per-run output dir; override via env var, else fresh timestamped folder

# ============================================================
# HELPER FUNCTIONS
# ============================================================
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

run_python() { python "$@" || { echo "ERROR: python $1 failed (exit $?)" >&2; exit 1; }; }
run_r()      { Rscript "$@" || { echo "ERROR: Rscript $1 failed (exit $?)" >&2; exit 1; }; }

# Upload a local file to Google Drive; non-fatal on failure
upload() {
    local src="$1" dest="$2"
    rclone -v copy "$src" "$GDRIVE_BASE/$dest" --drive-shared-with-me \
        || echo "WARNING: upload of $src failed" >&2
}

trap 'echo "ERROR: Pipeline failed at line $LINENO (exit code $?)" >&2' ERR

# ============================================================
# Working directories and run log
# ============================================================
work_dir=$(pwd)
scripts="$(pwd)/scripts"

# Tee stdout+stderr into a run log under the output dir. Guarded so that when a
# subscript is invoked by the wrapper (which has already set up the tee), we
# don't pile on a second tee that double-writes every line.
mkdir -p "${OUTPUT_DIR}"
RUN_LOG="${OUTPUT_DIR}/run.log"
if [[ "${GENDIV_LOG_SETUP:-0}" != "1" ]]; then
    exec > >(tee -a "$RUN_LOG") 2>&1
    export GENDIV_LOG_SETUP=1
fi
log "Output dir: ${OUTPUT_DIR}"
log "Pipeline log: ${RUN_LOG}"
