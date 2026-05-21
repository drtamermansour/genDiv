#!/usr/bin/env bash
# Install NeEstimator v2.x Linux CLI binary (Do et al. 2014; doi:10.1111/1755-0998.12157).
#
# NeEstimator estimates contemporary effective population size from a
# single SNP sample via the LD method (Waples 2006; Waples & Do 2008,
# 2010). McGivney et al. 2020 used this tool on global Thoroughbreds
# (Ne = 330); using it on our Standardbred cohort gives a method-matched
# comparison.
#
# Source: https://github.com/bunop/NeEstimator2.X
# Maintained until 2019; the canonical Molecular Fisheries Laboratory
# distribution is v2.1 (Dec 2017); the bunop fork's latest release tag
# is v2.0.1 of the same source.
#
# Default install location: ${REPO_ROOT}/tools/NeEstimator/Ne2L
# Override with: NEESTIMATOR_INSTALL_DIR=/path/to/dir bash install_neestimator.sh
set -eo pipefail

NEESTIMATOR_VERSION="${NEESTIMATOR_VERSION:-v2.0.1}"
NEESTIMATOR_REPO="${NEESTIMATOR_REPO:-https://github.com/bunop/NeEstimator2.X.git}"
NEESTIMATOR_INSTALL_DIR="${NEESTIMATOR_INSTALL_DIR:-$(cd "$(dirname "$0")/../.." && pwd)/tools/NeEstimator}"
NEESTIMATOR_SRC_DIR="$NEESTIMATOR_INSTALL_DIR/src"
NE2L_BIN="$NEESTIMATOR_INSTALL_DIR/Ne2x"

mkdir -p "$NEESTIMATOR_INSTALL_DIR"

# Note: we build from source rather than downloading the pre-built v2.0.1
# Ne2L Linux binary because the latter has been observed to crash with a
# "double free detected in tcache 2" on modern glibc (Ubuntu 22.04 +).
# The same source compiles cleanly to a working `Ne2x` binary that
# produces correct output.
if [[ -d "$NEESTIMATOR_SRC_DIR/.git" ]]; then
    echo "[install_neestimator] existing source at $NEESTIMATOR_SRC_DIR — fetching"
    git -C "$NEESTIMATOR_SRC_DIR" fetch --tags --quiet
else
    echo "[install_neestimator] cloning $NEESTIMATOR_REPO -> $NEESTIMATOR_SRC_DIR"
    git clone --quiet "$NEESTIMATOR_REPO" "$NEESTIMATOR_SRC_DIR"
fi
git -C "$NEESTIMATOR_SRC_DIR" checkout --quiet "$NEESTIMATOR_VERSION"

if [[ ! -x "$NE2L_BIN" ]]; then
    echo "[install_neestimator] compiling Ne2x via make (gcc required)"
    if ! command -v gcc >/dev/null 2>&1; then
        echo "[install_neestimator] ERROR: gcc not in PATH; cannot build NeEstimator." >&2
        exit 1
    fi
    (cd "$NEESTIMATOR_SRC_DIR" && make >/dev/null 2>&1)
    if [[ -x "$NEESTIMATOR_SRC_DIR/Ne2x" ]]; then
        cp "$NEESTIMATOR_SRC_DIR/Ne2x" "$NE2L_BIN"
    fi
fi

if [[ ! -x "$NE2L_BIN" ]]; then
    echo "[install_neestimator] ERROR: build did not produce $NE2L_BIN" >&2
    exit 1
fi

echo "[install_neestimator] binary verified at $NE2L_BIN"
echo "[install_neestimator] export NE2L_BIN=\"$NE2L_BIN\" to use from run_neestimator.sh"
