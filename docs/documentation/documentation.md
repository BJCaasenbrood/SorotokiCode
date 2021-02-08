---
layout: default
title: Documentation
nav_order: 3
has_children: true
permalink: docs/documentation
---

# Documentation
The Sorotoki toolkit includes a library of class-oriented functions for designing, modeling, and controlling a wide variety of soft robots. To mitigate complexity, the toolkit is designed to solve a multi-physics soft robotic problem by decomposing it into smaller sub-problems. 

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
