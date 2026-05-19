#!/usr/bin/env bash
# Install GONE2 (Santiago, Köpke & Caballero 2025; doi:10.1038/s41467-025-61378-w).
#
# GONE2 is a per-generation N_e trajectory estimator from a single SNP
# sample. Source: https://github.com/esrud/GONE2
#
# This script clones the repo, checks out the latest stable release tag,
# and verifies the binary is executable. Re-running is safe: if the
# install directory already exists, the script only updates and
# re-validates.
#
# Default install location: ${REPO_ROOT}/tools/GONE2/
# Override with: GONE2_INSTALL_DIR=/path/to/dir bash install_gone2.sh
#
# After install, run_gone2.sh consumes GONE2_BIN (the path to the
# GONE2 executable) — set via environment or via the wrapper's CLI.
set -eo pipefail

GONE2_VERSION="${GONE2_VERSION:-v1.0.2}"
GONE2_REPO="${GONE2_REPO:-https://github.com/esrud/GONE2.git}"
GONE2_INSTALL_DIR="${GONE2_INSTALL_DIR:-$(cd "$(dirname "$0")/../../tools" && pwd)/GONE2}"

mkdir -p "$(dirname "$GONE2_INSTALL_DIR")"

if [[ -d "$GONE2_INSTALL_DIR/.git" ]]; then
    echo "[install_gone2] existing checkout at $GONE2_INSTALL_DIR — fetching"
    git -C "$GONE2_INSTALL_DIR" fetch --tags --quiet
else
    echo "[install_gone2] cloning GONE2 from $GONE2_REPO -> $GONE2_INSTALL_DIR"
    git clone --quiet "$GONE2_REPO" "$GONE2_INSTALL_DIR"
fi

git -C "$GONE2_INSTALL_DIR" checkout --quiet "$GONE2_VERSION"
echo "[install_gone2] checked out $GONE2_VERSION"

# GONE2 ships as C++ source; build with make.
# Default MAXLOCI=2,000,000 and MAXIND=2,000 are sufficient for SNP-array
# studies (we use ~46K LD-pruned SNPs and <600 animals). If your panel
# exceeds these defaults, override via:
#   make MAXLOCI=10000000 MAXIND=5000 gone
GONE2_BIN="$GONE2_INSTALL_DIR/gone2"
if [[ ! -x "$GONE2_BIN" ]]; then
    echo "[install_gone2] compiling GONE2 (g++ required)"
    if ! command -v g++ >/dev/null 2>&1; then
        echo "[install_gone2] ERROR: g++ not in PATH; cannot build GONE2."
        exit 1
    fi
    (cd "$GONE2_INSTALL_DIR" && make gone)
fi

if [[ ! -x "$GONE2_BIN" ]]; then
    echo "[install_gone2] ERROR: build did not produce $GONE2_BIN"
    exit 1
fi

echo "[install_gone2] binary verified at $GONE2_BIN"
echo "[install_gone2] export GONE2_BIN=\"$GONE2_BIN\" to use from run_gone2.sh"
