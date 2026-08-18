---
name: steady-streaming-postprocess
description: Compute, diagnose, export, and plot Eulerian, ALE, or particle-based Lagrangian steady-streaming velocity from periodic FSI or rigid-wall HDF5/XDMF snapshots. Use for cycle-mean velocity fields, one-cycle particle drift, fluid-subdomain XDMF output, 1 mm axial slice integration, cranial/caudal/net streaming strength, area-normalized mean streaming velocity, FSI/RW comparisons, elastic-modulus sweeps, and parallel SLURM-friendly post-processing.
---

# Steady Streaming Postprocess

Use the bundled scripts instead of rebuilding this workflow ad hoc. Read
`references/steady_streaming_definitions.md` before choosing Eulerian, ALE, or Lagrangian averaging,
or before interpreting differences between them.

## Inspect Inputs

1. Inspect HDF5 groups and dataset shapes before running. Confirm mesh geometry, tetrahedral
   topology, cell tags, velocity snapshots, displacement snapshots, and snapshot indices.
2. Confirm that selected snapshots cover exactly one converged cardiac cycle in chronological
   order. Exclude a duplicated periodic endpoint.
3. Use all available uniformly spaced snapshots unless the output cadence is intentionally
   oversampled. Apply a stride only when it preserves one uniformly sampled cycle. The bundled
   scripts do not implement nonuniform-time quadrature or variable particle time steps.
4. Keep units consistent. These scripts assume mesh coordinates in mm, velocity in mm/s, area in
   mm^2, and strength in mm^3/s; CSV output also provides mL/s.
5. Use Python with NumPy, h5py, PyVista, SciPy, and Matplotlib. Use the legacy `dolfin` environment
   for `--reader fenics-function`; that reader reconstructs CG1 vector subfunctions, including
   velocity written as a subspace of a mixed velocity-pressure function.
6. For particle tracking, supply physical fluid velocity in the laboratory frame, not velocity
   relative to the ALE mesh.

## Choose the Mean

- Use `scripts/compute_eulerian_steady_streaming.py --frame ale` to reproduce direct averaging of
  velocity degrees of freedom at fixed ALE node identities. Label the result ALE cycle mean.
- Use the same script with `--frame eulerian` to interpolate each moving-domain snapshot onto the
  fixed physical grid from the reference snapshot before averaging. Check the reported valid-point
  fractions; this comparison is defined only on the common spatial support. Do not send an
  Eulerian field containing nonfinite boundary vertices to slice integration. Restrict/remesh it to
  the common support first.
- Use `scripts/compute_lagrangian_steady_streaming_particles.py` for the Parras/Sanchez construction:
  seed particles at fluid nodes for multiple starting phases, advect each set for one period,
  compute `[x(T)-x(0)]/T`, map each drift from trajectory-mean positions to a common mesh, and
  average over starting phases. Prefer Heun integration. Inspect valid-particle fractions and IDW
  remapping distances; refine temporal sampling when displacement per step is large.

Never average velocity magnitude when a signed vector mean is required. Never describe direct ALE
nodal averaging as a strict fixed-spatial Eulerian mean on a moving mesh.

## Compute A Field

For legacy FEniCS snapshots in one HDF5 file:

```bash
python scripts/compute_eulerian_steady_streaming.py \
  --mesh-h5 mesh.h5 --input snapshots.h5 --reader fenics-function \
  --indices 0:64 --velocity-path '/velocity/{index}' \
  --displacement-path '/displacement/{index}' --frame ale \
  --output-h5 cycle_mean.h5 --output-xdmf cycle_mean.xdmf
```

For nodal `(N,3)` HDF5 arrays, use `--reader h5-array`. The input filename itself may contain
`{index}` for one-file-per-snapshot checkpoints.

For Lagrangian drift:

```bash
python scripts/compute_lagrangian_steady_streaming_particles.py \
  --mesh-h5 mesh.h5 --input snapshots.h5 --reader fenics-function \
  --indices 0:64 --period-s 0.8 --start-count 20 \
  --velocity-path '/velocity/{index}' --displacement-path '/displacement/{index}' \
  --output-csv particles.csv --output-vtp particles.vtp \
  --field-h5 lagrangian_mean.h5 --field-xdmf lagrangian_mean.xdmf
```

Omit `--displacement-path` for a rigid-wall mesh. The output mesh is the first-snapshot fluid mesh.
Each phase-specific drift is associated with the particle trajectory's mean position and remapped
to that common mesh before phase averaging.

## Compute Slice Strength

Run `scripts/compute_strength_by_slices.py` on an Eulerian, ALE, or Lagrangian field HDF5:

```bash
python scripts/compute_strength_by_slices.py \
  --field-h5 cycle_mean.h5 --fluid-tag 258514 --slice-spacing-mm 1 \
  --case FSI --label FSI --output-csv fsi_strength.csv
```

The script cuts zero-thickness planes normal to z. It triangulates each cut and integrates the
positive and negative parts of the linearly interpolated axial velocity exactly, including
mixed-sign triangles. It writes:

- cranial strength: `integral(max(v_z,0), dA)`;
- caudal strength: `integral(min(v_z,0), dA)`, signed negative;
- net strength: cranial plus caudal;
- directional mean velocity: each strength divided by section area;
- distance from FM: `z_max - z_slice`.

Use `--axial-sign -1` only when the model's coordinate orientation makes negative axial velocity
cranial. Slices are explicitly normal to global z. By default, distance is `z_max-z_slice`; pass
`--fm-z-mm` when the FM origin is a known different z coordinate.

## Plot Profiles

Use standardized CSVs with `scripts/plot_steady_streaming_strength.py`. For an FSI/RW signed plot:

```bash
python scripts/plot_steady_streaming_strength.py \
  --profile FSI=fsi_strength.csv --profile RW=rw_strength.csv \
  --mode signed --quantity strength --color-by direction \
  --profile-style FSI=solid --profile-style RW=dashed \
  --output-stem figures/fsi_vs_rw_signed_strength
```

The plotter emits PNG, PDF, and SVG without a title. Defaults use Nimbus Roman/STIX, 22 pt labels
and ticks, 18 pt legends, red cranial curves, blue caudal curves, distance in cm, and mL/s.

## Parallel Runs

- Parallelize ALE accumulation over snapshots with `--snapshot-chunk K/N --partial-output partK.npz`,
  then pass every part to one `--combine-partials` invocation. Part indices must exactly cover the
  requested cycle.
- Parallelize Lagrangian phase averaging with a SLURM array over `--start-run K`, using the same
  `--start-count N` and a shared `--phase-output-dir`. Omit explicit output paths in worker jobs;
  the script derives collision-safe `lss_phase_K.*` names from `--start-run`.

```bash
#SBATCH --array=0-19
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
python scripts/compute_lagrangian_steady_streaming_particles.py \
  --mesh-h5 mesh.h5 --input snapshots.h5 --reader fenics-function \
  --indices 0:64 --period-s 0.8 --velocity-path '/velocity/{index}' \
  --displacement-path '/displacement/{index}' --start-count 20 \
  --start-run "${SLURM_ARRAY_TASK_ID}" --phase-output-dir phase_fields
```

- After all phase tasks finish, run:

```bash
python scripts/average_lagrangian_phase_fields.py \
  --phase-dir phase_fields --indices 0:20 \
  --output-h5 lagrangian_mean.h5 --output-xdmf lagrangian_mean.xdmf
```

  This verifies common-mesh fingerprints and finite values before averaging. Do not average phase
  files by copying their vectors onto the undeformed reference mesh. It also verifies source mesh,
  snapshot indices, field paths, period, integrator, particle-validity fraction, and IDW remapping
  diagnostics. Add explicit validity/distance thresholds when appropriate.
- Parallelize slice integration with `--slice-chunk K/N`, write one CSV per task, then concatenate
  rows once and sort by `distance_from_fm_mm`.
- Prefer SLURM arrays across independent cases. Do not let multiple tasks write the same output.

## Validate Outputs

Check metadata JSON files, snapshot coverage, field units, fluid cell count, valid-point/particle
fractions, and the expected near-cancellation `Q_cranial + Q_caudal` for a closed section. Treat
near-zero net strength as a mass-balance check, not evidence that recirculatory streaming is absent.
Slice integration fails deliberately when the selected field contains NaN or infinite values.
When the snapshot count is not divisible by `--start-count`, starting phases are approximately
uniform with floor-spaced offsets; inspect `start_offsets` in metadata.
