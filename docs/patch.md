---
layout: default
title: Patch notes
nav_order: 7
---

# Patch notes

## SOROTOKI - Alpha - 2.24.01 - Jan 24 - 2022
- Fixed broken installer. `vernum.m` file was missing on Repo. It has been
  replaced with `soropatch.m` which also includes the patch notes.

- **Sdf.m**:
  - Added a new function `Sdf.showcontour()`. It will show the contour of 2D
  signed distance functions. Currently not implemented for 3D Sdfs.

- **Open issue (2.24.01)**:    
  - @martijnschouten Missing DOI for citation, and long-term support/access.

## SOROTOKI - Alpha - 2.13.01 - Jan 13 - 2022
- Moved SOROTOKI from early alpha to alpha (prepping for official release).

- **Shapes.m**:
  - Significant update to the class Shapes. Shapes now requests a Fem class
  , a Matlab function_handle of functionals, or a evaluated matrix of
  shape functions (); to construct the class. Then, it serves as an input
  for the class Model (e.g., dynamic modeling of soft robots), or as an
  Inverse Kinematic solver for beam models.
  - Added a function `Shapes.jointEstimate(Fem)`. This will produce an
  optimal set of modal coefficient to reconstruct the beam model from
  FEM-driven data (accessed by Fem.Log).
  - Current `Shapes.jointEstimate()` exploits the fact that the shapes are
  orthonormal w.r.t. their integral over the spatial domain [0,L]. If shapes do not satisfy this conditions, inaccurate estimates can occur.
  - Added `Y = gsog_poly(X)` which ensure the columns of Y are mutually
  orthonormal derived from the matrix X. X must be full-rank. This
  function is simply the gramm-smith method for the innerproduct space
  int_C y_i y_j ds on the domain C := [0,L].

## SOROTOKI - Early Alpha - 1.12.06 - Dec 6 - 2021
- **Fem.m**:
  - Fixed some numerical issues with `Fem.simulate()` function. Stability is
  now further improved for larger timesteps.
  - Fixed some numerical issues with `Fem.Contact`. A wider range of
  TimeSteps result now in stable solutions. Still the timestep size
  should be taken with care. Also fixed the issues of induced
  (unstable) oscillations due to fast impact. Contact is now based on the
  initial Modulus of the hyper-elastic material,i.e., `Material.Emod()`.
  - Added time-base function_handle for the external pressures in dynamic
  finite-element simulations. They can added using:
   `Fem.AddConstraint('Pressure',id, @(t) sin(t))`;
  - Dynamic simulations now record the potential and kinetic energies.
  - Future implementation will have Load and Tendon-based dynamic forces
  - The potential energy of the external load is still missing...

## SOROTOKI - Early Alpha - 1.12.02 - Dec 2 - 2021
- **Fem.m**:
  - Added Fem.simulate() function to the Nonlinear Finite Element. The
  function is a standard Newmark-Beta dynamic solver routine for NLFEM.
  Alternatively, `fem.simulate()` can be used as an alternative for
  fem.solve() by setting the Fem.Material.Zeta large (overdamped).
  - Fixed the stability issue of Fem.Contact. Contact can now be added with
  `Fem.AddConstraint('Contact',@(x) sdf(x), [x,y])` (x and y move the SDF).
  for best stability, use `Fem.simulate()` with a `TimeStep < 1/75`.

- **Open issue (1.12.02)**:
  - Fem.simulate does not use dynamic external forces, rather
  the forces are scaled with a sigmoid function f_ext = f*sigmoid(t);
  This will be removed in the future with `fext = @(t) f* ....` - an
  explicit function of time that can be evaluated at time t.
  - Possible implementations: Displace, Load, Pressure, Gravity.


## SOROTOKI - Early Alpha - 1.02.0 - Oct 6 - 2021
- Added patch.md file to keep track of any changes to SOROTOKI.
- **Mesh.m**:
    - Fixed an issue where Mesh.showSDF was not producing plots.

## SOROTOKI - Early Alpha - 1.01.0
- Official beta release of SOROTOKI
