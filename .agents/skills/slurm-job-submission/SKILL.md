---
name: slurm-job-submission
description: Create, review, validate, submit, and troubleshoot Slurm batch jobs for Purdue RCAC clusters such as Bell, Negishi, and Gautschi. Use for .sh/.sub/.sbatch scripts, FEniCS or Python MPI jobs, shared-memory post-processing, hybrid MPI/OpenMP layouts, job arrays, account/partition/QOS selection, CPU binding, thread oversubscription, interactive-versus-batch performance differences, and job monitoring. Submit, cancel, or requeue a live job only when the user explicitly requests that action.
---

# Slurm Job Submission

Build resource-consistent jobs and make the requested allocation match the program's actual parallel model.

## Workflow

1. Determine the cluster, workload type, account, partition/QOS, walltime, memory, environment, working directory, command, and expected output location. Inspect an existing successful script from the same project when available.
2. Classify the workload as pure MPI, shared-memory, hybrid MPI/threaded, GPU, or job array. Read [references/job-patterns.md](references/job-patterns.md) for the matching resource equations and launcher rules.
3. For Bell, Negishi, or Gautschi, read [references/rcac-clusters.md](references/rcac-clusters.md). Verify live partitions, modules, and MPI plugins when possible; cached hardware values are not a substitute for the scheduler.
4. Copy the nearest template from `assets/` and replace every `<placeholder>`. Preserve project-specific module and MPI choices when they are known to work. Remove optional mail directives when notifications are not requested. For GPU jobs, build from current cluster-specific documentation instead of applying the standard CPU-node templates.
5. Run both checks before submission:

   ```bash
   bash -n run.sbatch
   python scripts/validate_slurm_script.py run.sbatch --cluster gautschi
   ```

   If available, use `sbatch --test-only run.sbatch` for a scheduler-side dry run. Do not create production directories merely to make a draft or test-only check pass; skip that check and report the prerequisite instead. Treat warnings as prompts for deliberate review, not automatic failures.
6. Submit with `sbatch run.sbatch` only when explicitly asked. Return the job ID and the exact script path. Never cancel, requeue, or alter a live job without explicit authorization.

## Resource rules

- Pure MPI: set `--cpus-per-task=1`, use one rank per requested core, and cap `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `NUMEXPR_NUM_THREADS` at 1.
- Shared-memory code: set `--ntasks=1`; set `--cpus-per-task` to the thread/worker count. Use `OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK` only for actual OpenMP/threaded work. Keep it at 1 for Python multiprocessing and pass `$SLURM_CPUS_PER_TASK` to the worker pool.
- Hybrid code: require `tasks_per_node * cpus_per_task <= physical_cores_per_node`; set `OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK`.
- Full-node CPU counts are Bell 128, Negishi 128, and Gautschi 192 for their standard CPU nodes. Do not apply those counts to GPU nodes, and verify current hardware before generating a new production script.
- Do not equate maximum rank count with maximum throughput. Benchmark representative core counts; nonlinear solvers, memory bandwidth, collectives, and parallel I/O often stop scaling before a full node.

## Script requirements

- Use `#!/bin/bash -l` when using `source activate`, `conda activate`, arrays, or `set -o pipefail`; `/bin/sh` is not reliably Bash.
- Put every `#SBATCH` directive before the first executable shell command. Shell variables do not expand inside directives; Slurm substitutions such as `%x`, `%j`, `%A`, and `%a` do.
- Request `--account`, `--partition`, `--nodes`, tasks, CPUs per task, walltime, job name, stdout, and stderr explicitly. Add memory, QOS, mail, GPUs, constraints, or arrays only when needed.
- Use `--ntasks-per-node` for predictable multi-node placement. Avoid an ambiguous total `-n` when rank placement matters.
- Use `module --force purge`, load the exact project modules, activate the environment, print `module list`, and `cd -- "$SLURM_SUBMIT_DIR"` or an intentional absolute working directory.
- Do not mix an `mpirun` executable from one MPI installation with `mpi4py`, PETSc, DOLFIN, or another library built against a different MPI. Mirror a proven project script or verify the MPI vendor first.
- For pure MPI, use either a site-tested `srun --mpi=<plugin> --cpu-bind=cores --kill-on-bad-exit=1` command or `mpirun -np "$SLURM_NTASKS" --bind-to core`. Do not assume a PMIx plugin name without checking `srun --mpi=list`.
- Let `#SBATCH --output` and `--error` capture logs unless a separate application log is genuinely useful. The parent directories for Slurm log paths must exist before `sbatch` starts the job.
- Write large results to scratch or depot as appropriate. Create result directories inside the job and use quoted absolute paths.
- Inspect the application's output interface before adding command-line arguments. Never stage and rewrite application source inside a batch script to work around a hard-coded output path. Report the mismatch and change the application only when the user authorizes that source edit. Treat concurrent jobs targeting one fixed output directory as unsafe.

## Diagnose performance

Compare interactive and batch runs only when node type, rank/thread layout, environment, input, and filesystem are equivalent. Record `hostname`, `module list`, `which python`, `which mpirun`, Slurm allocation variables, and all thread-limit variables. Check for rank/thread oversubscription before changing solver settings.

After completion, inspect `sacct`, RCAC `jobinfo` when available, stderr, maximum memory, elapsed time, and exit state. A faster interactive run does not by itself prove CPU contention; different node types, filesystems, MPI builds, or problem sizes can produce the same symptom.

## Templates

- `assets/pure-mpi-python.sbatch`: FEniCS, mpi4py, or other one-thread-per-rank jobs.
- `assets/shared-memory-python.sbatch`: serial, threaded, or Python worker-pool post-processing.
- `assets/hybrid-mpi-openmp.sbatch`: multiple MPI ranks with multiple cores per rank.
- `assets/snapshot-array.sbatch`: collision-safe independent chunks with bounded concurrency.

Write generated scripts only to the destination the user requested. Default to the user's project for real work and to a temporary directory for evaluations or dry runs; never place generated artifacts inside this skill.
