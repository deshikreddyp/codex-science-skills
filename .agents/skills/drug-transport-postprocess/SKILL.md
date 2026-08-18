---
name: drug-transport-postprocess
description: Post-process CSF drug-transport simulations with PyVista to compute cross-sectional net concentration c_net and area-average concentration c_avg versus time and caudal distance from the foramen magnum. Use for result*.xdmf folders or compact FEniCS concentration/displacement HDF5 snapshots, true-FSI warping or rigid-wall/CMM geometry, Paper2-compatible wide CSVs, long CSV/NPZ/metadata, heatmaps, time profiles, MPI extraction, and SLURM snapshot-chunk merging.
---

# Drug Transport Post-process

Extract Paper2-compatible `c_net` and `c_avg` on fixed axial planes while keeping simulation and depot inputs read-only. Treat `/depot/hgomezdi/data/dputluru/Paper2_figures/Transport` as the primary convention source.

## Workflow

1. Establish input type, requested snapshots/times, output directory, concentration units, and any supplied anatomical FM coordinate.
2. Read [input formats](references/input_formats.md) and [metric definitions](references/metric_definitions.md).
3. Inspect before processing:

   ```bash
   python "$SKILL_DIR/scripts/inspect_transport_inputs.py" --results-dir RESULTS
   python "$SKILL_DIR/scripts/inspect_transport_inputs.py" --mesh-h5 MESH.h5 --snapshots-h5 SNAPSHOTS.h5
   ```

   Confirm concentration, optional displacement/domain fields, topology, tags, indices, embedded times, coordinate scale, and concentration units.
4. Process and plot. Example:

   ```bash
   python "$SKILL_DIR/scripts/postprocess_concentration_slices.py" \
     --results-dir RESULTS --output-dir OUTPUT --profile-times-s 1800,3600
   ```

   For compact HDF5, use `--snapshots-h5` and `--mesh-h5`, select indices as needed, and supply `--dt-s` or `--period-s`. Override the inferred FM plane with `--fm-z-mm` when an anatomical coordinate is known.
5. Verify `net_concentration.csv`, `avg_concentration.csv`, long CSV, compressed NPZ, selected-timestep CSV, metadata JSON, and PNG/PDF plots. Inspect NaN patterns.

## Required choices

- Threshold fluid when labels exist. Warp by displacement for true FSI. Use stored rigid-wall/CMM geometry when displacement is absent and record the choice.
- Hold first-snapshot physical z planes fixed through time.
- Preserve integer distances `0, 1, 2, ... mm`, increasing caudally from the FM, with zero at the top of plots.
- Compute `c_net = integral_A c dA` and `c_avg = c_net/A`.
- Always label `c_net` as a cross-sectional concentration integral, not global 3D mass.
- Preserve concentration units and label `c_net` as concentration times mm2.
- Leave result/depot folders read-only. Refuse existing output deliverables unless `--overwrite` is explicit.

## Parallel processing

Use MPI over snapshots for large XDMF series. For compact HDF5 SLURM arrays, run serial `--snapshot-chunk K/N` jobs and merge only with:

```bash
python "$SKILL_DIR/scripts/merge_transport_partials.py" PARTIALS_ROOT --output-dir OUTPUT
```

Require the merge fingerprint, field-path, slice-coordinate, time, and value checks to pass.

## Handoff

Report snapshot/time selection, FM and slice extent, geometry warping decision, units, NaN slices, plot limits, and generated files.
