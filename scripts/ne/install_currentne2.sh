#!/usr/bin/env bash
# Install currentNe2 (Santiago, Köpke & Caballero 2025; doi:10.1038/s41467-025-61378-w).
#
# currentNe2 estimates contemporary effective population size from SNP
# data using LD between (mostly unlinked) loci. It is the single-point
# contemporary analogue to GONE2 and the recommended replacement for
# NeEstimator v2 in the same Nat Comms paper:
#   "GONE2, for inferring recent changes in Ne when a genetic map is
#    available, and currentNe2, which estimates contemporary Ne even in
#    the absence of genetic maps."
#
# Source: https://github.com/esrud/currentNe2 (C++ + OpenMP; same group as
# GONE/GONE2).
#
# Default install location: ${REPO_ROOT}/tools/currentNe2/
# Override with CURRENTNE2_INSTALL_DIR / CURRENTNE2_REPO env vars.
set -eo pipefail

CURRENTNE2_REPO="${CURRENTNE2_REPO:-https://github.com/esrud/currentNe2.git}"
CURRENTNE2_VERSION="${CURRENTNE2_VERSION:-main}"
CURRENTNE2_INSTALL_DIR="${CURRENTNE2_INSTALL_DIR:-$(cd "$(dirname "$0")/../../tools" && pwd)/currentNe2}"

mkdir -p "$(dirname "$CURRENTNE2_INSTALL_DIR")"

if [[ -d "$CURRENTNE2_INSTALL_DIR/.git" ]]; then
    echo "[install_currentne2] existing checkout at $CURRENTNE2_INSTALL_DIR — fetching"
    git -C "$CURRENTNE2_INSTALL_DIR" fetch --tags --quiet
else
    echo "[install_currentne2] cloning $CURRENTNE2_REPO -> $CURRENTNE2_INSTALL_DIR"
    git clone --quiet "$CURRENTNE2_REPO" "$CURRENTNE2_INSTALL_DIR"
fi
git -C "$CURRENTNE2_INSTALL_DIR" checkout --quiet "$CURRENTNE2_VERSION"

CURRENTNE2_BIN="$CURRENTNE2_INSTALL_DIR/currentne2"
if [[ ! -x "$CURRENTNE2_BIN" ]]; then
    echo "[install_currentne2] compiling currentNe2 (g++ + OpenMP required)"
    if ! command -v g++ >/dev/null 2>&1; then
        echo "[install_currentne2] ERROR: g++ not in PATH; cannot build currentNe2." >&2
        exit 1
    fi
    (cd "$CURRENTNE2_INSTALL_DIR" && make >/dev/null)
fi

if [[ ! -x "$CURRENTNE2_BIN" ]]; then
    echo "[install_currentne2] ERROR: build did not produce $CURRENTNE2_BIN" >&2
    exit 1
fi

echo "[install_currentne2] binary verified at $CURRENTNE2_BIN"
echo "[install_currentne2] export CURRENTNE2_BIN=\"$CURRENTNE2_BIN\" to use from run_currentne2.sh"
