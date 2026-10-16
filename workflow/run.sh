#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
/home/koki/anaconda3/envs/snakemake/bin/snakemake \
    --snakefile workflow/Snakefile \
    --cores "${1:-4}" \
    --printshellcmds \
    --keep-going \
    --rerun-incomplete \
    "${@:2}"
