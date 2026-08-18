# Steady-Streaming Definitions and Numerical Conventions

## Contents

1. Coordinate and sign convention
2. Eulerian, ALE, and Lagrangian means
3. Snapshot quadrature
4. Particle advection
5. Slice-wise streaming strength
6. Interpretation and checks

## 1. Coordinate and Sign Convention

Let the anatomical axial coordinate be `z`, with larger `z` toward the foramen magnum (FM). Define

\[
d = z_{\max}-z,
\]

so `d=0` is at the FM and distance increases caudally. With this convention, positive `v_z` is
cranial and negative `v_z` is caudal. If a mesh uses the opposite orientation, multiply the axial
component by `-1` before separating directions.

## 2. Eulerian, ALE, and Lagrangian Means

For a period `T`, the true Eulerian cycle mean at a fixed physical point is

\[
\overline{\mathbf u}_E(\mathbf x)
= \frac{1}{T}\int_{t_0}^{t_0+T}\mathbf u(\mathbf x,t)\,dt.
\]

On a deforming domain, this requires interpolation to fixed physical locations that remain in the
common spatial support of all snapshots.

Let the ALE map be `x=chi(X,t)`. Directly averaging a velocity degree of freedom attached to the
same mesh identity computes

\[
\overline{\mathbf u}_{ALE}(\mathbf X)
= \frac{1}{T}\int_{t_0}^{t_0+T}
\mathbf u(\boldsymbol\chi(\mathbf X,t),t)\,dt.
\]

This is an ALE pullback mean. It is not exactly the fixed-spatial Eulerian mean, although the two
can be close for small mesh motion.

For a particle flow map `Phi_T`, a one-cycle Lagrangian drift velocity from phase `t_i` is

\[
\overline{\mathbf u}_L(\mathbf x_0)
= \frac{\boldsymbol\Phi_T(\mathbf x_0)-\mathbf x_0}{T}.
\]

It contains the net stroboscopic particle drift, including the interaction between oscillatory
velocity and spatial gradients. A phase-independent Parras/Sanchez-style field repeats this
calculation for several `t_i`, associates each drift with its trajectory-mean position, maps each
result to a common Eulerian mesh, and averages over `t_i`. In perturbation language, Lagrangian mean
drift differs from the Eulerian mean by Stokes drift after both quantities are represented in a
compatible frame.

## 3. Snapshot Quadrature

For `N` uniformly spaced periodic snapshots with no duplicated endpoint,

\[
\overline{\mathbf u}\approx \frac{1}{N}\sum_{n=0}^{N-1}\mathbf u_n.
\]

Do not include both `t=0` and `t=T` when they contain the same phase. For nonuniform times, use
time-interval weights rather than an arithmetic mean; the bundled scripts require uniform spacing
and must be adapted for nonuniform times. A stride is valid only when it leaves a uniform
representation of one complete cycle. For example, 640 snapshots may be reduced with a stride only
when the selected phases still cover the full period uniformly; a stride of 3 is not a generic
correction for having roughly three times as many snapshots.

## 4. Particle Advection

Choose uniformly distributed starting phases over the buffered cycle. For each phase, seed
particles at the fluid vertices of that phase's physical mesh. At each time, construct the current
physical mesh using displacement and interpolate physical fluid velocity at particle positions.
Use laboratory-frame fluid velocity, not fluid velocity relative to the ALE mesh.
For equal step size `dt=T/N`, Heun's method uses

\[
\mathbf x^*=\mathbf x_n+\Delta t\,\mathbf u(\mathbf x_n,t_n),
\]

\[
\mathbf x_{n+1}=\mathbf x_n+\frac{\Delta t}{2}
\left[\mathbf u(\mathbf x_n,t_n)+\mathbf u(\mathbf x^*,t_{n+1})\right].
\]

Wrap the last step to the selected starting phase. For every valid trajectory, also compute

\[
\mathbf x_o=\frac{1}{N+1}\sum_{n=0}^{N}\mathbf x_n.
\]

Treat `(x_o, [x(T)-x(0)]/T)` as scattered vector data, interpolate it to the common output mesh by
local inverse-distance weighting, and average the mapped fields over starting phases. This is the
workflow used in `/home/dputluru/fenics_run/steady_streaming_better_toy`.

For complex SAS geometry, inspect the nearest-source remapping distances. Euclidean IDW can mix
opposite sides of a thin wall if the source cloud is too sparse or displaced; reduce the neighbor
count, refine temporal/particle sampling, or use topology-aware remapping when this occurs.

Mark a trajectory invalid if interpolation leaves the fluid mesh. Boundary seeds can be
numerically delicate; report the valid fraction and nearest-source remapping distances. A single
starting phase is a stroboscopic drift diagnostic, not the final phase-averaged LSS field.

## 5. Slice-Wise Streaming Strength

For an axial section `Gamma_d`, define signed directional strength from any steady-streaming field:

\[
Q_{SS}^{cra}(d)=\int_{\Gamma_d}\max(v_{SS,z},0)\,dA,
\]

\[
Q_{SS}^{cau}(d)=\int_{\Gamma_d}\min(v_{SS,z},0)\,dA.
\]

The caudal value is stored as a negative number. The net is

\[
Q_{SS}^{net}=Q_{SS}^{cra}+Q_{SS}^{cau}.
\]

Area-normalized directional velocities are `Q/A`. They describe section-wide exchange strength,
not the local peak streaming velocity.

Use a zero-thickness planar cut. Triangulate the cut and integrate the positive and negative parts
of the piecewise-linear nodal field. On a mixed-sign triangle, locate the zero crossings and
integrate each sign region; assigning the entire triangle according to its centroid or mean creates
a mesh-dependent bias.

## 6. Interpretation and Checks

- Incompressible recirculatory flow in a closed section should have approximately
  `Q_cra + Q_cau = 0`, apart from geometric, interpolation, and solver error.
- Equal and opposite flow strengths do not imply zero solute flux because concentration need not
  be equal in the two streaming branches.
- Lagrangian and Eulerian/ALE means should not be expected to match. Their difference is physically
  meaningful and can be sensitive to phase, spatial gradients, vortices, and wall motion.
- A noisy mean usually indicates wrong snapshot ordering, an incomplete cycle, a duplicated phase,
  incompatible degree-of-freedom ordering, or averaging fields without first placing them in the
  intended common frame.
- When phase workers are used, verify identical geometry, topology, and cell-tag fingerprints
  before averaging their mapped Lagrangian fields.
- Always save the selected indices, period, reader mode, field paths, fluid tag, slice spacing, and
  sign convention beside outputs.
