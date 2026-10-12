#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
snakemake \
    --snakefile workflow/Snakefile \
    --cores "${1:-4}" \
    --use-conda \
    --printshellcmds \
    --keep-going \
    --rerun-incomplete \
    "${@:2}"
