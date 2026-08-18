# Drug-transport metric definitions and conventions

The primary specification is `/depot/hgomezdi/data/dputluru/Paper2_figures/Transport`.

## Axial distance convention

- Use fixed global-z planes with 1 mm spacing by default.
- Preserve the established integer-distance convention: distances are `0, 1, 2, ... mm` increasing caudally from the FM.
- By default, set the FM plane to `floor(z_max)` of the first selected physical fluid geometry and the last plane to `ceil(z_min)`. Override with `--fm-z-mm` and `--z-bottom-mm` when the anatomical FM coordinate is known.
- Determine the planes from the first selected physical geometry once and hold them fixed through time.

## Metrics

- `c_net(d,t) = integral_A c dA`.
- `c_avg(d,t) = c_net(d,t)/A(d,t)`.

`c_net` is a cross-sectional concentration integral with units of concentration times area. It is not global 3D drug mass. Empty/nonintersecting slices are NaN rather than zero.

Triangulate each plane slice and integrate a CG1 point field exactly as triangle area times the mean of the three vertex values.

## Outputs and plots

Preserve `net_concentration.csv` and `avg_concentration.csv`: the first column is `time_seconds`, remaining header entries are times, and each data row begins with distance from the FM. Also write long-form CSV, compressed NPZ, selected-timestep CSV, and metadata JSON.

Default Paper2-style heatmaps put time on x and distance from FM on y, with zero at the top. Produce PNG and PDF. Use a warm sequential colormap and Nimbus Roman when available. Requested times produce distance profiles using the nearest sampled physical time.
