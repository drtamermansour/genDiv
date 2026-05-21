#!/usr/bin/env bash
# Install SNeP v1.11 (Barbato et al. 2015; doi:10.3389/fgene.2015.00109).
#
# SNeP estimates a per-generation N_e trajectory from LD decay using the
# Sved-Feldman (1971) approximation
#     Ne_t = (1 / (4 c_t)) * (1 / E[r^2]_t - alpha)
# where c_t maps to generations ago t via the chosen mapping function
# (Haldane by default), E[r^2]_t is the mean pairwise r^2 in the c-bin,
# and alpha is a sample-size correction. SNeP is the dominant LD-decay
# N_e tool in the 2023-2025 livestock literature (see Manunza 2025).
#
# Source: https://sourceforge.net/projects/snepnetrends/  (static Linux
# x86_64 binary, no compilation required). The Linux build is named
# SNeP_111 (Ubuntu 20.x, GCC >= 4.8.2).
#
# Default install location: ${REPO_ROOT}/tools/SNeP/SNeP1.11
# Override with SNEP_INSTALL_DIR / SNEP_DOWNLOAD_URL env vars.
set -eo pipefail

SNEP_INSTALL_DIR="${SNEP_INSTALL_DIR:-$(cd "$(dirname "$0")/../../tools" && pwd)/SNeP}"
SNEP_BIN_NAME="${SNEP_BIN_NAME:-SNeP1.11}"
SNEP_BIN="$SNEP_INSTALL_DIR/$SNEP_BIN_NAME"
SNEP_DOWNLOAD_URL="${SNEP_DOWNLOAD_URL:-https://sourceforge.net/projects/snepnetrends/files/binaries/SNeP_111/download}"

mkdir -p "$SNEP_INSTALL_DIR"

if [[ -x "$SNEP_BIN" ]]; then
    echo "[install_snep] existing binary at $SNEP_BIN — already installed"
    exit 0
fi

# SNeP is distributed as a single statically-linked Linux binary on
# the snepnetrends SourceForge project (no archive). SourceForge's
# /files/<path>/download URL returns an HTML mirror-picker page with a
# signed `ts=<token>` redirect URL embedded as a meta-refresh. We scrape
# that signed URL and follow it with the picker page set as Referer.
echo "[install_snep] downloading SNeP from $SNEP_DOWNLOAD_URL"
if ! command -v curl >/dev/null 2>&1; then
    echo "[install_snep] ERROR: 'curl' not in PATH; cannot fetch from SourceForge" >&2
    exit 1
fi

picker_html=$(curl -sL -A "Mozilla/5.0" "$SNEP_DOWNLOAD_URL") || picker_html=""
signed_url=$(printf "%s" "$picker_html" \
              | grep -oE 'url=https://downloads\.sourceforge\.net/project/snepnetrends/binaries/SNeP_111\?[^"]+' \
              | head -1 \
              | sed -e 's/^url=//' -e 's/&amp;/\&/g')

if [[ -z "$signed_url" ]]; then
    echo "[install_snep] ERROR: could not find signed mirror URL in the picker page." >&2
    echo "[install_snep] Manually download SNeP_111 from $SNEP_DOWNLOAD_URL" >&2
    echo "[install_snep] and place it at: $SNEP_BIN" >&2
    exit 1
fi

if ! curl -sL -A "Mozilla/5.0" -e "$SNEP_DOWNLOAD_URL" -o "$SNEP_BIN" "$signed_url"; then
    echo "[install_snep] ERROR: download failed via signed mirror URL." >&2
    rm -f "$SNEP_BIN"
    exit 1
fi

# SourceForge sometimes serves an HTML interstitial instead of the file
# when the mirror picker fails. ELF binaries start with the 7F 45 4C 46
# magic ("\x7fELF"); reject anything else and tell the user.
if ! head -c 4 "$SNEP_BIN" | od -An -c | grep -q '177   E   L   F'; then
    head -c 200 "$SNEP_BIN" | tr -d '\0' > "$SNEP_BIN.headpeek" 2>/dev/null || true
    echo "[install_snep] ERROR: downloaded file is not an ELF binary (likely an HTML interstitial)." >&2
    echo "[install_snep] First bytes: $(head -c 60 "$SNEP_BIN" | tr -d '\0')" >&2
    echo "[install_snep] Manually download SNeP_111 from $SNEP_DOWNLOAD_URL" >&2
    echo "[install_snep] and place it at: $SNEP_BIN" >&2
    rm -f "$SNEP_BIN" "$SNEP_BIN.headpeek"
    exit 1
fi

chmod +x "$SNEP_BIN"

if [[ ! -x "$SNEP_BIN" ]]; then
    echo "[install_snep] ERROR: $SNEP_BIN is not executable after install" >&2
    exit 1
fi

echo "[install_snep] binary installed at $SNEP_BIN"
echo "[install_snep] export SNEP_BIN=\"$SNEP_BIN\" to use from run_snep.sh"
