---
name: csf-flow-postprocess
description: Post-process CSF FSI, rigid-wall, or CMM simulation results with PyVista to compute axial flow rate, cross-sectional area and relative area deformation, area-weighted pressure, and maximum tissue displacement. Use for result*.xdmf series or compact FEniCS mesh/snapshot HDF5 inputs, fixed z-slice heatmaps, requested-distance waveforms, requested-time spatial profiles, MPI extraction, and SLURM snapshot-chunk merging. Do not use for steady streaming, cycle-mean velocity, or particle drift; use steady-streaming-postprocess instead.
---

# CSF Flow Post-process

Compute reproducible CSF slice metrics without modifying simulation inputs. Treat `/depot/hgomezdi/data/dputluru/Paper2_figures/Flow` as the primary convention source.

## Workflow

1. Identify XDMF-series or compact FEniCS-HDF5 input. Establish requested metrics, snapshot selection, output directory, physical-time mapping for HDF5, and any nondefault units or tags.
2. Read [input formats](references/input_formats.md) and [metric definitions](references/metric_definitions.md).
3. Inspect before processing:

   ```bash
   python "$SKILL_DIR/scripts/inspect_csf_inputs.py" --results-dir RESULTS
   python "$SKILL_DIR/scripts/inspect_csf_inputs.py" --mesh-h5 MESH.h5 --snapshots-h5 SNAPSHOTS.h5
   ```

   Confirm fields, topology, domain tags, time range, snapshot indices, and assumed units. Resolve mismatches through the documented overrides.
4. Run only the requested metrics. Example:

   ```bash
   python "$SKILL_DIR/scripts/postprocess_csf_slices.py" \
     --results-dir RESULTS --output-dir OUTPUT \
     --metrics flow-rate area-deformation pressure max-displacement \
     --waveform-distances-mm 10.5,40.5 --profile-times-s 0.5,1.0
   ```

   For compact HDF5, use `--snapshots-h5` and `--mesh-h5`, select indices if needed, and supply `--dt-s` or `--period-s`.
5. Verify the long CSV, compact NPZ, selected-snapshot CSV, metadata JSON, PNG files, and PDF files. Inspect NaN patterns and unavailable-metric metadata.

## Required choices

- Threshold fluid before slicing. Warp fluid and solid by displacement for true FSI. Use stored geometry when displacement is absent and record that choice.
- Use fixed global-z planes at 1 mm spacing, with the first center 0.5 mm below the first selected physical fluid `z_top`.
- Use the first selected snapshot as the area-deformation baseline.
- Preserve pressure units unless the user supplies both an explicit scale and unit label.
- Treat a single-domain mesh without a domain field as fluid. Never infer tissue displacement when no tissue domain exists.
- Leave result and depot folders read-only. Refuse existing output deliverables unless `--overwrite` is explicit.

## Parallel processing

Use MPI over snapshots for large XDMF series. For compact HDF5 SLURM arrays, run serial jobs with `--snapshot-chunk K/N`, then merge only with:

```bash
python "$SKILL_DIR/scripts/merge_csf_partials.py" PARTIALS_ROOT --output-dir OUTPUT
```

Require the merge fingerprint, field-path, slice-coordinate, time, and value checks to pass.

## Handoff

Report input files, selected times, slice range, warped/static geometry, units and scales, unavailable metrics, NaN slices, and generated files.

Redirect cycle-mean Eulerian/ALE/Lagrangian streaming, streaming strength, and particle-drift requests to `$steady-streaming-postprocess` unchanged.
