#!/bin/bash -l
#SBATCH --job-name=<job_name>
#SBATCH --account=<account>
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<workers>
#SBATCH --time=<HH:MM:SS>
#SBATCH --mem=<memory>
#SBATCH --output=<existing_absolute_log_dir>/%x.o%j
#SBATCH --error=<existing_absolute_log_dir>/%x.e%j
#SBATCH --mail-user=<email>
#SBATCH --mail-type=FAIL,END
#SBATCH --export=ALL

set -euo pipefail

module --force purge
module load <anaconda_module>
set +u
source activate <conda_environment>
set -u
module list

# This template assumes a Python process pool. Keep each worker's BLAS/OpenMP pool at one.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

cd -- "$SLURM_SUBMIT_DIR"
mkdir -p "<absolute_results_dir>"

python3 -u "<absolute_script.py>" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --output-dir "<absolute_results_dir>"
