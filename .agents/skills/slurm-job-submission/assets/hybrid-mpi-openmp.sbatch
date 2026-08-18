#!/bin/bash -l
#SBATCH --job-name=<job_name>
#SBATCH --account=<account>
#SBATCH --partition=cpu
#SBATCH --nodes=<nodes>
#SBATCH --ntasks-per-node=<mpi_ranks_per_node>
#SBATCH --cpus-per-task=<threads_per_rank>
#SBATCH --time=<days-HH:MM:SS>
#SBATCH --output=<existing_absolute_log_dir>/%x.o%j
#SBATCH --error=<existing_absolute_log_dir>/%x.e%j
#SBATCH --mail-user=<email>
#SBATCH --mail-type=FAIL,END
#SBATCH --export=ALL

set -euo pipefail

module --force purge
module load <compiler_and_mpi_modules>
module list

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export NUMEXPR_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cd -- "$SLURM_SUBMIT_DIR"

# Replace this with a launcher verified for the loaded MPI. It must reserve
# SLURM_CPUS_PER_TASK processing elements per rank; pure-MPI binding is insufficient.
<verified_hybrid_launcher> "<absolute_executable>" <arguments>
