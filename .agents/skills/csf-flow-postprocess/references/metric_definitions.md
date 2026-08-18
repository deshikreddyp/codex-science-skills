# CSF metric definitions and conventions

The primary specification is the archived Paper2 Flow output under
`/depot/hgomezdi/data/dputluru/Paper2_figures/Flow`. Preserve these conventions unless the user explicitly requests a different definition.

## Fixed axial slices

- Work in global coordinates with planes normal to z.
- Default spacing is 1 mm.
- The first slice center is 0.5 mm below the top of the first selected physical fluid geometry.
- `distance_from_FM = z_top(first selected physical fluid geometry) - z_plane`.
- Determine the planes once and reuse them for every snapshot. Do not move the planes with the tissue.
- Coordinates and displacement are assumed to be mm. Velocity is assumed to be mm/s.

## Metrics

- `flow_rate_mm3_per_s = integral_A v_z dA`. Positive values point toward increasing global z.
- `area_mm2 = integral_A 1 dA` on the physical fluid slice.
- `relative_area_deformation = (A(d,t)-A(d,t0))/A(d,t0)`, where `t0` is the first selected snapshot, not necessarily simulation time zero.
- `mean_pressure = integral_A p dA / A`. Preserve solver pressure units unless `--pressure-scale` and `--pressure-unit` are supplied.
- `max_tissue_displacement_mm = max_A_solid |u|` on the solid-domain plane slice.

Do not infer a tissue maximum when a solid domain or its displacement is unavailable. Record it as unavailable. A rigid-wall/CMM dataset without displacement uses stored geometry and has zero relative area deformation for a time-invariant mesh.

## Integration

Triangulate each plane slice and integrate a CG1 point field exactly per triangle using triangle area times the arithmetic mean of its three vertex values. Empty/nonintersecting slices are NaN, not zero.

## Plotting

Default heatmaps use time on x and distance from FM on y, with zero distance at the top. Store flow in mm3/s but label the default heatmap in mL/s after multiplying by 1e-3. Use Nimbus Roman if available, with a serif fallback. Record any explicit color limits or unit scaling in metadata.
