# CSF input formats and reliability rules

## Inspect first

Run `scripts/inspect_csf_inputs.py` before extraction. Confirm fields, centers, topology, domain values, physical times, snapshot indices, coordinate scale, and supplied solver units. Never modify a result or depot folder.

## XDMF series

Expected layout is `result*.xdmf` plus referenced HDF5 data. Defaults are point fields `Velocity`, `Pressure`, and `Displacement`, and cell field `Sub-domain`. XDMF visualization output normally uses fluid tag 0 and solid tag 1. Use embedded XDMF physical times; fail if they are absent rather than inventing times. The computation script distributes snapshots across MPI ranks.

## Compact FEniCS HDF5

Supply `--mesh-h5` and `--snapshots-h5`. Defaults are `/velocity/INDEX`, `/pressure/INDEX`, and `/displacement/INDEX`, with native fluid tag 258514 and solid tag 258515. Supply `--dt-s` or `--period-s` when physical time metadata is unavailable.

The reader reconstructs FEniCS functions with `dolfin.HDF5File.read()` and maps degrees of freedom to vertices with `vertex_to_dof_map`. Never treat raw `vector_0` storage order as nodal ordering.

Use `--snapshot-chunk K/N` for serial SLURM arrays. Each chunk writes under `partial_K_of_N`; merge only with `merge_csf_partials.py`, which verifies source fingerprint, field paths, slice convention/coordinates, unique times, and the absence of infinities. NaNs are permitted only as explicit nonintersecting/unavailable slice values. Use `--overwrite` deliberately; otherwise existing deliverables are protected.

All field names, groups, domain tags, FM coordinate, spacing, index selection, time spacing/period, and output location have CLI overrides. For a single-domain rigid-wall mesh with no domain array, treat the full stored mesh as fluid; do not invent a solid domain.
