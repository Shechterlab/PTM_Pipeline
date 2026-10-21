#!/bin/bash

set -eo pipefail

CODE_ROOT=${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}

if command -v rg >/dev/null 2>&1; then
  SEARCH_CMD=(rg -n '^(<<<<<<<|=======|>>>>>>>)' "$CODE_ROOT")
else
  SEARCH_CMD=(grep -RInE '^(<<<<<<<|=======|>>>>>>>)' "$CODE_ROOT")
fi

# grep and rg exit 1 when nothing matches and 2 on an error: only 1 passes
markers_status=0
MARKERS=$("${SEARCH_CMD[@]}" 2>&1) || markers_status=$?
if [ "$markers_status" -eq 0 ]; then
  echo "[ptm_hygiene] merge conflict markers detected:"
  printf '%s\n' "$MARKERS"
  exit 1
elif [ "$markers_status" -ne 1 ]; then
  echo "[ptm_hygiene] the conflict-marker search failed:"
  printf '%s\n' "$MARKERS"
  exit 1
fi

if git -C "$CODE_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  # Frozen inputs and release files are tracked on purpose: each is listed, with its checksum, in a manifest.
  MANIFESTED=$(
    [ -f "$CODE_ROOT/inputs/MANIFEST.tsv" ] && tail -n +2 "$CODE_ROOT/inputs/MANIFEST.tsv" | cut -f1 | sed 's|^|inputs/|'
    for sums in "$CODE_ROOT"/releases/*/MD5SUMS; do
      [ -f "$sums" ] || continue
      dir=$(dirname "${sums#"$CODE_ROOT"/}")
      awk '{print $2}' "$sums" | sed "s|^|$dir/|"
    done)
  FORBIDDEN_TRACKED=$(git -C "$CODE_ROOT" ls-files \
    '.codex' '.codex/*' '.agents/*' 'scripts.zip' '*.zip' '*.pyc' 'scripts/__pycache__/*' \
    'PTM_data/*' 2>/dev/null | while IFS= read -r file; do
      status=$(git -C "$CODE_ROOT" status --porcelain -- "$file" || true)
      if printf '%s\n' "$status" | grep -Eq '^( D|D )'; then
        continue
      fi
      if printf '%s\n' "$MANIFESTED" | grep -Fxq "$file"; then
        continue
      fi
      if [ -e "$CODE_ROOT/$file" ]; then
        printf '%s\n' "$file"
      fi
    done || true)
  if [ -n "$FORBIDDEN_TRACKED" ]; then
    echo "[ptm_hygiene] generated/local artifacts are tracked and should be removed from Git:"
    printf '%s\n' "$FORBIDDEN_TRACKED"
    exit 1
  fi

  LARGE_TRACKED=$(git -C "$CODE_ROOT" ls-files -z | while IFS= read -r -d '' file; do
    path="$CODE_ROOT/$file"
    if printf '%s\n' "$MANIFESTED" | grep -Fxq "$file"; then
      continue
    fi
    if [ -f "$path" ]; then
      size=$(wc -c < "$path")
      if [ "$size" -gt 10000000 ]; then
        printf '%s\t%s bytes\n' "$file" "$size"
      fi
    fi
  done)
  if [ -n "$LARGE_TRACKED" ]; then
    echo "[ptm_hygiene] tracked files larger than 10 MB detected:"
    printf '%s\n' "$LARGE_TRACKED"
    exit 1
  fi
fi

echo "[ptm_hygiene] repo hygiene check passed"
