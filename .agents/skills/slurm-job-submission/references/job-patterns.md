# Slurm job patterns

## Pure MPI

Use for DOLFIN/FEniCS, PETSc, mpi4py, and compiled MPI programs that intentionally distribute work across ranks.

```text
cpus_per_rank = 1
ranks_per_node <= physical_cores_per_node
total_ranks = nodes * ranks_per_node
```

Cap numerical-library thread pools at one thread. Pin ranks to cores. More ranks can increase Krylov communication and parallel-I/O overhead, so benchmark rather than assuming a full node is fastest.

Use one launcher whose MPI ABI matches the application:

```bash
srun --mpi=<verified_plugin> --ntasks="$SLURM_NTASKS" --cpu-bind=cores --kill-on-bad-exit=1 python3 -u solver.py
# or, when the environment is known to support OpenMPI's launcher:
mpirun -np "$SLURM_NTASKS" --bind-to core python3 -u solver.py
```

## Shared-memory or worker-pool Python

Use one Slurm task and multiple CPUs per task:

```text
ntasks = 1
cpus_per_task = desired threads or workers
```

For OpenMP/BLAS threading, set `OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK` and normally set the other library caps to the same value. For Python multiprocessing or `ProcessPoolExecutor`, keep BLAS/OpenMP caps at one and pass `$SLURM_CPUS_PER_TASK` as the worker count; otherwise every worker can create its own thread pool.

## Hybrid MPI plus threads

Use only when the application is designed and tested for hybrid execution:

```text
tasks_per_node * cpus_per_task <= physical_cores_per_node
OMP_NUM_THREADS = cpus_per_task
```

MPI binding must reserve `cpus_per_task` processing elements for each rank. Launcher syntax varies across OpenMPI, MVAPICH2, and Slurm plugins; copy a tested hybrid launcher instead of applying pure-MPI `--bind-to core` blindly.

## Job arrays

Use arrays for independent simulations or snapshot chunks. Set a concurrency cap with `%`, such as `--array=0-15%4`, to limit simultaneous I/O and scheduler load. Give each task a unique directory based on `$SLURM_ARRAY_TASK_ID`; never let array tasks overwrite the same partial output.

After all chunks succeed, submit a separate merge job with an `afterok` dependency:

```bash
array_id=$(sbatch --parsable chunks.sbatch)
sbatch --dependency="afterok:${array_id}" merge.sbatch
```

Do not submit this chain unless the user explicitly requests live submission.

## Output and reproducibility

- Slurm opens `#SBATCH --output` and `--error` before the shell runs, so their parent directories must already exist at submission time.
- Use `%x` for job name, `%j` for job ID, and `%A_%a` for array job/task IDs.
- Print `date`, `hostname`, `module list`, the executable paths, `SLURM_JOB_ID`, `SLURM_NTASKS`, `SLURM_CPUS_PER_TASK`, and relevant thread variables at job start.
- Use unbuffered Python (`python3 -u` or `PYTHONUNBUFFERED=1`) for useful live logs.
- Do not write all MPI ranks to separate uncoordinated files unless the application is designed for it.

## Common failure modes

- **Batch is much slower than interactive:** first check node type, task/thread counts, MPI vendor, CPU affinity, filesystem, and input equality. Oversubscription is common but not the only cause.
- **Job stays pending:** inspect `squeue -j <jobid> -o '%.18i %.9P %.8T %.10M %.6D %R'` and `scontrol show job <jobid>` for the reason.
- **Immediate environment failure:** compare `module list`, `which python`, and package/MPI paths with the working interactive environment.
- **Out of memory:** inspect `sacct -j <jobid> --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode` or RCAC `jobinfo <jobid>`; request memory deliberately rather than a full node by habit.
- **Poor scaling:** compare elapsed time and core-hours over a small strong-scaling sweep. Optimize time-to-solution and allocation efficiency, not CPU occupancy alone.

Useful same-node starting sweeps are 32/64/128 ranks on Bell or Negishi and 48/96/192 ranks on Gautschi. Keep the node type, environment, input, output filesystem, and solver tolerances identical. Request a realistic walltime for each test because unnecessarily long limits can increase queue wait.
