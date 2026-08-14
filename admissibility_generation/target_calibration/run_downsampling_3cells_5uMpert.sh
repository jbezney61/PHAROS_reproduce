#!/bin/bash
#SBATCH --job-name=rsc_python
#SBATCH --output=logs/rsc_python_%j.out
#SBATCH --error=logs/rsc_python_%j.err

set -euo pipefail

cd /oak/stanford/groups/larsms/Users/jbezney/tahoe100m

export APPTAINER_CACHEDIR=/oak/stanford/groups/larsms/Users/jbezney/Software/singularity_cache
export APPTAINER_TMPDIR=/oak/stanford/groups/larsms/Users/jbezney/Software/singularity_tmp
export SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR
export SINGULARITY_TMPDIR=$APPTAINER_TMPDIR

export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16

module load apptainer

mkdir -p logs

apptainer exec \
  --bind $PWD:/workspace \
  rsc.sif \
  bash -lc "cd /workspace && python downsample_tahoe_3cells_5um_plus_dmso_cpu.py"

