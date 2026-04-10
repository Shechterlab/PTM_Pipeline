#!/bin/bash

CODE_ROOT=${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}

if rg -n '^(<<<<<<<|=======|>>>>>>>)' "$CODE_ROOT" >/tmp/ptm_conflict_markers.txt 2>/dev/null; then
  echo "[ptm_hygiene] merge conflict markers detected:"
  cat /tmp/ptm_conflict_markers.txt
  exit 1
fi

echo "[ptm_hygiene] repo hygiene check passed"
