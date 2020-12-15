---
layout: default
title: Documentation
nav_order: 3
has_children: true
permalink: docs/documentation
---

# Documentation
The Sorotoki toolkit includes a library of class-oriented functions for designing, modeling, and controlling a wide variety of soft robots. To mitigate complexity, the toolkit is designed to solve a multi-physics soft robotic problem by decomposing it into smaller sub-problems. 
<!-- For instance, the classes [`Mesh.m`](./documentation/meshing.html) and [`Fem.m`](./finite-elements.html) provide some numerical tools to model the continuum material behavior of soft structures, whose extrapolated material models can later be used in the dynamical modeling with the class [`Dynamic.m`](./documentation/dynamic-modeling.html).  -->

## Class-oriented functions
The main classes of Sorotoki are listed below:

| Class   | Matlab command     | Description  |
| *Mesh.lib* | [`msh = Mesh()`](./documentation/meshing.html) | Meshing class |
| *Fem.lib* | [`fem = Fem()`](./documentation/meshing.html) | (Nonlinear) finite elements class |
| *Magnet.lib* | [`mfem = Mfem()`](./documentation/meshing.html) | Magnetics finite elements class |
| *Model.lib* | [`mdl = Model()`](./documentation/dynamic-modeling.html) | Dynamical modeling class |
| *Gmodel.lib* | [`obj = Gmodel()`](./documentation/graphics.html) | Computer-graphical model class |
| *Rig.lib* | [`rig = Rig()`](./documentation/graphics.html) | Inverse-Kinematics (IK) rigging class |
| *Control.lib* | [`ctr = Control()`](./documentation/graphics.html) | (Real-time) control class |
{: .fs-3}

<!-- ```matlab
% list of classes
msh = Mesh();	    % meshing class
fem = Fem();   	    % finite elements class
mag = Mfem();       % magnetics finite elements class
obj = Gmodel();     % graphical model class
mdl = Model();      % dynamical model class
rig = Rig();        % IK rigging class
ctr = Control();    % (real-time) control class
``` -->
