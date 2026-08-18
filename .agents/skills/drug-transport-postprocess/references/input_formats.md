# Drug-transport input formats and reliability rules

Run `scripts/inspect_transport_inputs.py` first. Confirm topology, `Concentration`, optional `Displacement`, optional `Sub-domain`, domain tags, physical times, coordinate scale, concentration unit, and selected indices. Existing simulation and depot folders remain read-only.

## XDMF result folders

Read `result*.xdmf` and referenced HDF5 data. Defaults are point scalar `Concentration`, point vector `Displacement`, cell field `Sub-domain`, and XDMF fluid tag 0. Threshold to fluid when a domain field is present. Warp by displacement when present; otherwise use stored rigid-wall/CMM geometry and record that choice. Require embedded physical times. XDMF snapshots can be distributed with MPI.

## Compact FEniCS HDF5

Supply `--mesh-h5` and `--snapshots-h5`. Defaults are `/concentration/INDEX`, optional `/displacement/INDEX`, and native fluid tag 258514. Supply `--dt-s` or `--period-s`.

Reconstruct every field using `dolfin.HDF5File.read()` and `vertex_to_dof_map`. Raw `vector_0` ordering is never nodal ordering.

For SLURM arrays, use `--snapshot-chunk K/N`; each chunk writes `partial_K_of_N`. Merge with `merge_transport_partials.py`, which checks source fingerprint, field paths, slice convention/coordinates, unique times, and infinities. NaN nonintersecting slices are allowed. Existing outputs are refused unless `--overwrite` is explicit.

Field names/groups, fluid tag, FM and bottom coordinates, spacing, snapshots, time mapping, units, plot limits, profile times, and output directory are all overrideable.
