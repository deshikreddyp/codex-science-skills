# Purdue RCAC cluster reference

Treat this as a cached starting point. Before creating a production script, verify the active cluster, partitions, QOS limits, modules, and node type against the scheduler and current RCAC documentation.

## Standard CPU nodes

| Cluster | Physical CPU cores/node | Typical memory/node | CPU architecture | Recommended MPI family |
|---|---:|---:|---|---|
| Bell | 128 | 256 GB | 2 x 64-core AMD EPYC Rome | OpenMPI |
| Negishi | 128 | 256 GB | 2 x 64-core AMD EPYC Milan | OpenMPI or MVAPICH2 |
| Gautschi | 192 | 384 GB | 2 x 96-core AMD EPYC Genoa | OpenMPI |

These numbers describe standard CPU nodes. High-memory nodes can have more memory, and GPU nodes can have different CPU counts. For example, Gautschi's documented GPU nodes have 112 or 128 CPU cores rather than 192. Bell is documented with a 2026 retirement date, so verify that it is still available before targeting it.

## Live checks

Use short, read-only checks; do not run compute work on a login node.

```bash
hostname -f
slist
mybalance
timeout 15s sinfo -h -o '%P|cores=%c|memory_MB=%m|gres=%G' | sort -u
srun --mpi=list
module spider anaconda
module spider openmpi
```

- Use the account shown by `slist`/`mybalance`, not a guessed account.
- The usual CPU partition is `cpu`; request `--partition=cpu` explicitly.
- Use `--qos=standby` only when the user accepts its current walltime and priority/preemption policy. Query the live cluster documentation because those limits differ by cluster and can change.
- Use `module --force purge` followed by the same module versions used to build the application. The user's recent Gautschi FEniCS jobs use an Anaconda module plus `source activate fenicsproject`; keep that proven pattern unless the environment is intentionally changed.
- Check an MPI-aware Python environment with `which python`, `which mpirun`, and `python -c 'from mpi4py import MPI; print(MPI.get_vendor())'` before choosing a launcher from a different module tree.

## Official sources

Verified 2026-08-18:

- Bell overview: https://db.rcac.purdue.edu/knowledge/bell/overview
- Bell running jobs: https://docs.rcac.purdue.edu/userguides/bell/run_jobs/
- Negishi overview: https://rcac.purdue.edu/compute/negishi
- Gautschi overview: https://db.rcac.purdue.edu/compute/gautschi
- Gautschi accounts, partitions, and QOS: https://docs.rcac.purdue.edu/userguides/gautschi/run_jobs/queues/
- RCAC Conda jobs: https://docs.rcac.purdue.edu/userguides/gilbreth/run_jobs/examples/python_conda/

Prefer the cluster-specific equivalent of a general RCAC page when one exists.
