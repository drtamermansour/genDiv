#!/usr/bin/env bash
set -eo pipefail

# Load CONFIG, helper functions, working-dir vars, and log redirection.
source "$(dirname "$0")/genDiversity_common.sh"


## Run whole-pop preprocessing (Sections 1–3).
bash "$(dirname "$0")/genDiversity_shared.sh"

## Per-group stage. Each invocation produces reference files for one group.
for rg in wholePop Trotter Pacer; do
    bash "$(dirname "$0")/genDiversity_per_group.sh" "$rg"
done

############################################
## x. Nucleotide diversity statistic (pi) -- This section is under development
############################################
## Nucleotide diversity is a population-level metric, the average number of differences between a pair of chromosomes, across all chromosome combinations within the population.
## This is distinct from simply measuring heterozygosity.
## Variant vs. Invariant Sites: Traditional pi calculations require knowledge of both variant and invariant sites (i.e., sequencing data).
##    With a SNP array, pi will be overestimated because it ignores the conserved (non-variable) parts of the genome (i.e., it is "SNP-based" diversity rather than a true "genomic" diversity.)


## Cross-group aggregation.
bash "$(dirname "$0")/genDiversity_aggregate.sh"
