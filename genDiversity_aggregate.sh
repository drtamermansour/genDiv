#!/usr/bin/env bash
# Cross-group aggregation stage for the genDiv pipeline.
# Runs once after genDiversity_per_group.sh has been invoked for all three
# groups. Contains only the blocks that truly require all three per-group
# runs: twoGait / threeBooksize per-sample concatenations, and the
# Froh-vs-ROHsh plots that consume them. Per-group wholePop-only analyses
# (F_ROH histograms, roh_high, gait/bookSize stratified summaries) live in
# genDiversity_per_group.sh under the if [[ "$rg" == "wholePop" ]] blocks.
set -eo pipefail
source "$(dirname "$0")/genDiversity_common.sh"

log "Cross-group aggregation stage"

roh_RG="${OUTPUT_DIR}/divStats/roh.L3"

##########################################
## FST statistics (needs all 3 per-group iterations done)
##########################################
## fst_stats.R reads all five FST summary files — three whole-pop ones
## produced in per_group.sh wholePop §3a (sex / gait / bookSize) plus two
## per-gait book-size ones produced in per_group.sh §3 for Trotter / Pacer.
Rscript scripts/fst_stats.R "${OUTPUT_DIR}" &> ${OUTPUT_DIR}/divStats/fst_stats.txt
rclone -v copy ${OUTPUT_DIR}/divStats/fst_stats.txt "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Fst/" --drive-shared-with-me

##########################################
## Cross-group per-sample ROH_sh concatenations
##########################################
## twoGait   = Trotter + Pacer rows from each gait's own consensus-intersect file.
## threeBooksize = the six book-size subgroup rows (Trotter_LOW/MEDIUM/HIGH +
##                 Pacer_LOW/MEDIUM/HIGH).
## Each sample therefore appears exactly once in each concatenation, under the
## consensus of the gait (or book-size subgroup) it actually belongs to.
head -n1 "${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt" \
    > "${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt"
for rg in Trotter Pacer; do
    tail -n+2 "${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt"
done >> "${roh_RG}.perSample_intersect_twoGait_consensus_${pct}pct.summary.txt"

head -n1 "${roh_RG}.perSample_intersect_wholePop_consensus_${pct}pct.summary.txt" \
    > "${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt"
for rg in Trotter_LOW Trotter_MEDIUM Trotter_HIGH Pacer_LOW Pacer_MEDIUM Pacer_HIGH; do
    tail -n+2 "${roh_RG}.perSample_intersect_${rg}_consensus_${pct}pct.summary.txt"
done >> "${roh_RG}.perSample_intersect_threeBooksize_consensus_${pct}pct.summary.txt"

##########################################
## F_ROH vs ROH_sh plots (consume twoGait/threeBooksize summaries built above)
##########################################
## One plot per aggregation view. $froh is the gait_bookSize-stratified F_ROH
## table produced by per_group.sh wholePop (§8b); by the time aggregate.sh
## runs, that file already exists.
froh="${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt"
for gp in wholePop twoGait threeBooksize; do
    conShare="${roh_RG}.perSample_intersect_${gp}_consensus_${pct}pct.summary.txt"
    output_file="${OUTPUT_DIR}/divStats/Froh_vs_ROHsh_${gp}.png"
    python scripts/roh_plot.py "$conShare" "$froh" "$output_file"
    rclone -v copy "$output_file" --drive-shared-with-me "remote_UCDavis_GoogleDr:STR_Imputation_2025/outputs/Froh/"
    output_prefix="${OUTPUT_DIR}/divStats/normalized_ROHsh_${gp}"
    run_python scripts/roh_histograms.py --metric ratio "$conShare" "$froh" "$output_prefix"
    upload "$output_prefix".histogram.png "Froh/"
    upload "$output_prefix".density.png   "Froh/"
    output_prefix2="${OUTPUT_DIR}/divStats/ROHshared_${gp}"
    run_python scripts/roh_histograms.py --metric shared "$conShare" "$froh" "$output_prefix2"
    upload "$output_prefix2".histogram.png "Froh/"
done &> "${OUTPUT_DIR}/divStats/roh_sh.log"

##########################################
## LD-decay effective population size (Ne)
##########################################
## Modular: each Ne tool wrapper lives at scripts/ne/run_<tool>.sh and is
## invoked iff its wrapper exists and is executable. Any individual tool
## can be removed by deleting (or chmod -x) its wrapper file; no other
## edits are required. Outputs (per-tool subdir + a standardised
## Ne_<tool>_summary.csv) land under ${OUTPUT_DIR}/divStats/ne/<tool>/.
##
## GONE2 (Santiago et al. 2025) consumes the post-QC but NOT LD-pruned
## PLINK set because its estimator uses the full LD spectrum across
## recombination-distance bins; LD pruning would discard the signal.
## NeEstimator and SNeP (when wired) use the LD-pruned set, matching
## McGivney 2020 / Manunza 2025 standard practice.
ne_dir="${OUTPUT_DIR}/divStats/ne"
mkdir -p "$ne_dir"
filtered_unpruned_prefix="${OUTPUT_DIR}/filtered/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered"
filtered_pruned_prefix="${OUTPUT_DIR}/LD_pruned/USTA_Diversity_Study.remap.refAlleles.dedup.plink1.filtered.norm.phased.LD_prune"

for ne_tool in gone2 currentne2 neestimator snep; do
    wrapper="$(dirname "$0")/scripts/ne/run_${ne_tool}.sh"
    [[ -x "$wrapper" ]] || continue
    case "$ne_tool" in
        gone2)
            bash "$wrapper" \
                --unpruned-prefix "$filtered_unpruned_prefix" \
                --group "wholePop:${OUTPUT_DIR}/preprocess/samples.wholePop.txt" \
                --group "Trotter:${OUTPUT_DIR}/preprocess/samples.Trotter.txt" \
                --group "Pacer:${OUTPUT_DIR}/preprocess/samples.Pacer.txt" \
                --out-dir "$ne_dir" \
                || log "WARNING: Ne tool ${ne_tool} failed (continuing)"
            ;;
        currentne2|neestimator|snep)
            bash "$wrapper" \
                --pruned-prefix "$filtered_pruned_prefix" \
                --group "wholePop:${OUTPUT_DIR}/preprocess/samples.wholePop.txt" \
                --group "Trotter:${OUTPUT_DIR}/preprocess/samples.Trotter.txt" \
                --group "Pacer:${OUTPUT_DIR}/preprocess/samples.Pacer.txt" \
                --out-dir "$ne_dir" \
                || log "WARNING: Ne tool ${ne_tool} failed (continuing)"
            ;;
    esac
    # Upload outputs for this tool: combined summary CSV at the top level,
    # plus per-group raw trajectories, summary STATS, and the recombination-
    # bin LD diagnostic + tool log for downstream debugging / reviewer
    # verification.
    if [[ -d "$ne_dir/${ne_tool}" ]]; then
        for f in "$ne_dir/${ne_tool}"/*.csv; do
            [[ -f "$f" ]] && upload "$f" "Ne/${ne_tool}/"
        done
        # Tool-specific small diagnostic files:
        #   GONE2:       *_GONE2_Ne   *_GONE2_STATS   *_GONE2_d2
        #   NeEstimator: input.*Ne.txt  input.*NexLD.txt  info.*.txt  options.*.txt
        #   SNeP:        *.NeAll   *.LDAll   *SNeP.log
        #   currentNe2:  *_currentNe2_OUTPUT.txt
        #   All tools:   *.log
        for f in "$ne_dir/${ne_tool}"/*/*_Ne \
                 "$ne_dir/${ne_tool}"/*/*_STATS \
                 "$ne_dir/${ne_tool}"/*/*_d2 \
                 "$ne_dir/${ne_tool}"/*/*Ne.txt \
                 "$ne_dir/${ne_tool}"/*/*NexLD.txt \
                 "$ne_dir/${ne_tool}"/*/info.*.txt \
                 "$ne_dir/${ne_tool}"/*/options.*.txt \
                 "$ne_dir/${ne_tool}"/*/*.NeAll \
                 "$ne_dir/${ne_tool}"/*/*.LDAll \
                 "$ne_dir/${ne_tool}"/*/*_currentNe2_OUTPUT.txt \
                 "$ne_dir/${ne_tool}"/*/*.log; do
            [[ -f "$f" ]] && upload "$f" "Ne/${ne_tool}/$(basename "$(dirname "$f")")/"
        done
    fi
done

##########################################
## ROH_common cross-group concatenations and Figures 2-4
##########################################
## Each sample is scored against its own group's landscape during per_group.sh.
## Aggregate.sh stitches those per-group TSVs into the twoGait and
## threeBooksize views (every sample appears exactly once), then drives
## the three ROH_common figures.
roh_common_dir="${OUTPUT_DIR}/divStats/roh_common"

# Per-sample concatenations (header from wholePop; rows from the per-gait /
# per-booksize files). Each sample appears exactly once in each view.
head -n1 "${roh_common_dir}/roh_common.wholePop.tsv" \
    > "${roh_common_dir}/roh_common.twoGait.tsv"
for rg in Trotter Pacer; do
    tail -n+2 "${roh_common_dir}/roh_common.${rg}.tsv"
done >> "${roh_common_dir}/roh_common.twoGait.tsv"

# threeBooksize: per_group.sh emits scoring TSVs only for the main gait,
# not for the six book-size subs. We build the threeBooksize view by
# joining the gait-level scores with the gait_bookSize factor file so
# each row carries its book-size label, exactly as the existing
# Froh_vs_ROHsh aggregation already does for ROH_share.
awk 'BEGIN{FS=OFS="\t"} \
     NR==FNR { bs[$2]=$3; next } \
     FNR==1 { print $0, "gait_bookSize" } \
     FNR>1  { print $0, (($1 in bs) ? bs[$1] : "undefined") }' \
    "${OUTPUT_DIR}/preprocess/USTA_Diversity_Study.gait_bookSize" \
    "${roh_common_dir}/roh_common.twoGait.tsv" \
    > "${roh_common_dir}/roh_common.threeBooksize.tsv"

upload "${roh_common_dir}/roh_common.twoGait.tsv"       "ROH/roh_common/"
upload "${roh_common_dir}/roh_common.threeBooksize.tsv" "ROH/roh_common/"

# Pairwise Mann-Whitney + Cohen's d across the six gait x book_size
# subgroups for genome-wide ROH_common + 4 length-class scores.
pairwise_stats="${roh_common_dir}/roh_common_pairwise_stats.tsv"
run_python scripts/roh_common_pairwise_stats.py \
    --roh-common "${roh_common_dir}/roh_common.twoGait.tsv" \
    --froh       "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt" \
    --out        "$pairwise_stats"
upload "$pairwise_stats" "ROH/roh_common/"

# Per-subgroup mean +/- SD summary (9 rows: wholePop + 2 gaits + 6 gait x book-size).
# wholePop uses cohort-wide landscape scoring; other rows use within-gait scoring.
subgroup_summary="${roh_common_dir}/roh_common_subgroup_summary.csv"
run_python scripts/roh_common_subgroup_summary.py \
    --wholepop      "${roh_common_dir}/roh_common.wholePop.tsv" \
    --threebooksize "${roh_common_dir}/roh_common.threeBooksize.tsv" \
    --out           "$subgroup_summary"
upload "$subgroup_summary" "ROH/roh_common/"

{
    # Figure 2 — raw-line Manhattan landscape per group, all five panels
    # (genome-wide + 4 length classes) in one PNG.
    for rg in wholePop Trotter Pacer; do
        out="${roh_common_dir}/Fig2_manhattan.${rg}.png"
        run_python scripts/roh_common_plot.py manhattan \
            --landscape-dir "$roh_common_dir" --rg "$rg" --out "$out"
        upload "$out" "ROH/roh_common/"
    done

    # Figure 3 — F_ROH vs ROH_common scatter (Pacer / Trotter panels).
    fig3="${roh_common_dir}/Fig3_FROH_vs_ROHcommon.png"
    run_python scripts/roh_common_plot.py scatter \
        --roh-common "${roh_common_dir}/roh_common.twoGait.tsv" \
        --froh       "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt" \
        --out        "$fig3"
    upload "$fig3" "ROH/roh_common/"

    # Figure 4 — length-stratified boxplots by book-size.
    fig4="${roh_common_dir}/Fig4_lengthbox.png"
    fig4b="${roh_common_dir}/Fig4b_lengthbox.png"
    run_python scripts/roh_common_plot.py lengthbox \
        --roh-common "${roh_common_dir}/roh_common.twoGait.tsv" \
        --froh       "${OUTPUT_DIR}/divStats/roh.L3_Froh_gait_bookSize.txt" \
        --out        "$fig4" \
        --nested-out "$fig4b"
    upload "$fig4"  "ROH/roh_common/"
    upload "$fig4b" "ROH/roh_common/"
} &> "${OUTPUT_DIR}/divStats/roh_common.log"

##########################################
## ROH islands: selection-signature candidates from the ROH_common landscapes
##########################################
## For each group (wholePop, Trotter, Pacer), threshold the per-window
## f_w landscape at the top-1% (with absolute f_w >= 0.5 flagged as
## high-confidence), merge contiguous high-f_w windows allowing small
## gaps, drop islands narrower than 500 kb, and annotate each island
## with overlapping Ensembl EquCab3 protein-coding genes. A curated
## horse-selection candidate-gene list flags islands that overlap
## known selection loci (DMRT3, MSTN, LCORL/NCAPG, MC1R, KIT, ASIP,
## STX17, MITF, ...). A cross-group consolidation pass identifies
## shared (multi-group) vs gait-specific islands.
roh_islands_dir="${OUTPUT_DIR}/divStats/roh_islands"
mkdir -p "$roh_islands_dir"
gtf="$(dirname "$0")/input_data/annotation/Equus_caballus.EquCab3.0.115.gtf.gz"
candidates="$(dirname "$0")/scripts/horse_selection_candidates.tsv"
if [[ -f "$gtf" && -f "$candidates" ]]; then
    run_python scripts/roh_islands_annotate.py \
        --group   "wholePop:${roh_common_dir}/landscape.wholePop.tsv" \
        --group   "Trotter:${roh_common_dir}/landscape.Trotter.tsv" \
        --group   "Pacer:${roh_common_dir}/landscape.Pacer.tsv" \
        --gtf     "$gtf" \
        --candidates "$candidates" \
        --top-percentile     1.0 \
        --absolute-threshold 0.5 \
        --max-gap-windows    2 \
        --min-island-kb      500 \
        --differential-pair  "Pacer:Trotter" \
        --differential-top-percentile 1.0 \
        --out-dir "$roh_islands_dir"
    for f in "$roh_islands_dir"/roh_islands.*.csv; do
        upload "$f" "ROH/roh_islands/"
    done
else
    log "WARNING: ROH islands step skipped — annotation files not found at"
    log "  $gtf"
    log "  $candidates"
fi
