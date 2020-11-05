---
layout: default
title: Installation and getting started
nav_order: 2
---

# Installation
{: .no_toc }

Download the latest stable version ([.zip](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master) or [.tar](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master)), and unpack the compressed folder at any desired directory. Alternatively, you can directly clone the repository with the command:
```git
git clone https://github.com/BJCaasenbrood/SorotokiCode.git
```

Setting up the toolkit is relatively straightforward. To set-up the toolkit, simply run the command below. That's it, the soft robotics toolkit is now ready-to-use. 
```rust
% installation command
sorotoki();
```
The command above is also used to update the toolkit. It is recommended to run `sorotoki.m` to check for updates occasionally. 

# Getting started
To get started, type the following line in the command window:
```rust
% show demos
sorotoki('demo');
```

The Sorotoki toolkit includes a library of class-oriented functions for designing, modeling, and controlling a wide variety of soft robots. To mitigate complexity, the toolkit is designed to solve a multi-physics soft robotic problem by decomposing it into smaller sub-problems. For instance, the classes [`Mesh`](./libary/meshing.html) and [`Fem`](./finite-elements.html) provide some numerical tools to model the continuum material behavior of soft structures, whose extrapolated material models can later be used in the dynamic modeling part with the class [`Mesh`](./libary/dynamic-modeling.html). 

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
