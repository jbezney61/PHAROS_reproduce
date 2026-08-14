#!/bin/bash
#SBATCH --job-name=preprocess_tahoe_embed
#SBATCH --output=logs/preprocess_tahoe_embed_%j.out
#SBATCH --error=logs/preprocess_tahoe_embed_%j.err

set -euo pipefail

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

export APPTAINER_CACHEDIR=/oak/stanford/groups/larsms/Users/jbezney/Software/singularity_cache
export APPTAINER_TMPDIR=/oak/stanford/groups/larsms/Users/jbezney/Software/singularity_tmp
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"

export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16

module load apptainer

apptainer exec \
  --bind "$PWD":/workspace \
  rsc.sif \
  bash -lc "cd /workspace && python preprocess_merged_tahoe.py \
    --input data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_raw_cpu.h5ad \
    --output data/merged_5um_perturbations_plus_DMSO_100_per_cell_line_log1p_norm10k.h5ad \
    --target-sum 10000" 