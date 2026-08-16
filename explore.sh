#!/usr/bin/env bash
# explore.sh — driver for one-off exploratory analyses on a finished pipeline run.
#
# Each task below is a self-contained block that calls a script in explore/.
# Tasks read from an existing $OUTPUT_DIR (a results_<timestamp> folder) and
# write any new artifacts back into it (typically under divStats/ or a
# task-specific subfolder). Nothing here is part of the main pipeline.
#
# Usage:
#   OUTPUT_DIR=results_20260421_200003 bash explore.sh [task_name ...]
#
# With no args, runs all tasks in order. With args, runs only the named tasks.
#
# Convention for new tasks:
#   - Outputs land in $OUTPUT_DIR/explore/ (local).
#   - Uploaded to "$GDRIVE_BASE/explore/<task_name>/" on Google Drive via
#     the `upload` helper (defined in genDiversity_common.sh). Every task
#     must rclone-upload its deliverables so they survive outside the
#     compute node.
#   - Each task block below carries a header listing: Question, Inputs,
#     Outputs (local), Remote (GDrive destination).

set -eo pipefail
source "$(dirname "$0")/genDiversity_common.sh"
module load rclone  ## matches genDiversity_shared.sh — needed for upload helper

tasks=("$@")
run_task() { [[ ${#tasks[@]} -eq 0 ]] || printf '%s\n' "${tasks[@]}" | grep -Fxq "$1"; }

##########################################
## relationship_comparison_colored
##########################################
## Question:
##   Does the Standard-GRM vs ROH-GRM pair-kinship relationship differ when
##   broken down by gait pairing or by book-size pairing, and does the pattern
##   shift between the 1Mb / 5Mb / 10Mb ROH windows?
##
## What it does:
##   Replicates the top-right "Relationship Comparison" subplot from
##   rep_ROHRM/roh_<cutoff>Mb.Threshold_3SD/Robust_Matrix_Comparison_Enhanced.wholePop.png
##   (scatter of Kinship_Std on x, Kinship_ROH on y, plus a 1:1 dashed line),
##   recolouring each pair by a categorical pair-label:
##     - Gait:      3 unordered categories (Pacer-Pacer, Pacer-Trotter,
##                  Trotter-Trotter). Pairs with an 'undefined' gait on
##                  either side are dropped.
##     - Book size: 6 unordered categories (HIGH-HIGH, HIGH-LOW, HIGH-MEDIUM,
##                  LOW-LOW, LOW-MEDIUM, MEDIUM-MEDIUM). Pairs where either
##                  sample is missing a book-size assignment are dropped.
##   Each colour scheme is a separate 1x3 figure, one panel per ROH cutoff,
##   axes shared across panels within a figure for direct visual comparison.
##   Legend entries include per-category pair counts.
##
##   The book-size view is emitted three times: once with all 6 categories
##   (full overview), and then two focused subsets that share the full-view
##   palette but re-fit their axes to the subset for better resolution of
##   fine structure:
##     - HIGHgroups:      HIGH-HIGH, HIGH-MEDIUM, HIGH-LOW
##     - LOWandMEDIUMgroups: LOW-LOW, LOW-MEDIUM, MEDIUM-MEDIUM
##
## Inputs (read from $OUTPUT_DIR):
##   rep_ROHRM/roh_1Mb.Threshold_3SD/Pairwise_Differences.wholePop.csv
##   rep_ROHRM/roh_5Mb.Threshold_3SD/Pairwise_Differences.wholePop.csv
##   rep_ROHRM/roh_10Mb.Threshold_3SD/Pairwise_Differences.wholePop.csv
##   preprocess/USTA_Diversity_Study.bookSize
##
## Outputs (local, under $OUTPUT_DIR/explore/):
##   relationship_comparison_byGait.wholePop.png                            (1x3, coloured by gait pair)
##   relationship_comparison_byBookSize.wholePop.png                        (1x3, all 6 book-size pair categories)
##   relationship_comparison_byBookSize_HIGHgroups.wholePop.png             (1x3, HIGH-* subset)
##   relationship_comparison_byBookSize_LOWandMEDIUMgroups.wholePop.png     (1x3, LOW/MEDIUM-* subset)
##
## Remote (GDrive, under $GDRIVE_BASE):
##   explore/relationship_comparison_colored/<all four pngs above>
if run_task relationship_comparison_colored; then
    log "Explore: relationship_comparison_colored"
    run_python explore/relationship_comparison_colored.py --output-dir "$OUTPUT_DIR"
    upload "${OUTPUT_DIR}/explore/relationship_comparison_byGait.wholePop.png"                         "explore/relationship_comparison_colored/"
    upload "${OUTPUT_DIR}/explore/relationship_comparison_byBookSize.wholePop.png"                     "explore/relationship_comparison_colored/"
    upload "${OUTPUT_DIR}/explore/relationship_comparison_byBookSize_HIGHgroups.wholePop.png"          "explore/relationship_comparison_colored/"
    upload "${OUTPUT_DIR}/explore/relationship_comparison_byBookSize_LOWandMEDIUMgroups.wholePop.png"  "explore/relationship_comparison_colored/"
fi
