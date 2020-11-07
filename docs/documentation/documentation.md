---
layout: default
title: Documentation
nav_order: 3
has_children: true
permalink: docs/documentation
---

# Documentation
The Sorotoki toolkit includes a library of class-oriented functions for designing, modeling, and controlling a wide variety of soft robots. To mitigate complexity, the toolkit is designed to solve a multi-physics soft robotic problem by decomposing it into smaller sub-problems. For instance, the classes [`Mesh.m`](./documentation/meshing.html) and [`Fem.m`](./finite-elements.html) provide some numerical tools to model the continuum material behavior of soft structures, whose extrapolated material models can later be used in the dynamic modeling part with the class [`Mesh.m`](./documentation/dynamic-modeling.html). 

The main classes of Sorotoki are listed below:

```matlab
% list of classes
msh = Mesh();	    % meshing class
fem = Fem();   	    % finite elements class
mag = Mfem();       % magnetics finite elements class
obj = Gmodel();     % graphical model class
mdl = Model();      % dynamical model class
rig = Rig();        % IK rigging class
ctr = Control();    % (real-time) control class
```
