#!/bin/bash -l
#SBATCH --job-name=<job_name>
#SBATCH --account=<account>
#SBATCH --partition=cpu
#SBATCH --nodes=<nodes>
#SBATCH --ntasks-per-node=<mpi_ranks_per_node>
#SBATCH --cpus-per-task=1
#SBATCH --time=<days-HH:MM:SS>
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

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

cd -- "$SLURM_SUBMIT_DIR"
mkdir -p "<absolute_results_dir>"

date --iso-8601=seconds
hostname -f
printf 'SLURM_JOB_ID=%s\nSLURM_NTASKS=%s\nSLURM_CPUS_PER_TASK=%s\n' \
  "$SLURM_JOB_ID" "$SLURM_NTASKS" "$SLURM_CPUS_PER_TASK"
command -v python3
command -v srun

# Verify the plugin with `srun --mpi=list`; use the value proven for this environment.
srun --mpi=<verified_pmix_plugin> --ntasks="$SLURM_NTASKS" --cpu-bind=cores --kill-on-bad-exit=1 \
  python3 -u "<absolute_script.py>" --output-dir "<absolute_results_dir>"
