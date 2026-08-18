#!/bin/bash -l
#SBATCH --job-name=<job_name>
#SBATCH --account=<account>
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=<first>-<last>%<max_concurrent>
#SBATCH --time=<HH:MM:SS>
#SBATCH --mem=<memory>
#SBATCH --output=<existing_absolute_log_dir>/%x_%A_%a.out
#SBATCH --error=<existing_absolute_log_dir>/%x_%A_%a.err
#SBATCH --mail-user=<email>
#SBATCH --mail-type=FAIL,END
#SBATCH --export=ALL

set -euo pipefail

module --force purge
module load <anaconda_module>
set +u
source activate <conda_environment>
set -u

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

cd -- "$SLURM_SUBMIT_DIR"
chunk_dir="<absolute_partial_root>/chunk_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$chunk_dir"

python3 -u "<absolute_script.py>" \
  --chunk-index "$SLURM_ARRAY_TASK_ID" \
  --chunk-count <total_chunks> \
  --output-dir "$chunk_dir"
