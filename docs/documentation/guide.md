---
layout: default
title: User guide
nav_order: 3
parent: Documentation
---

# User guide
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---

## Useful commands within SOROTOKI
<div class="code-example" markdown="1">
## `sorotoki(varargin)`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Function`{: .text-red-000} -- Calls the SOROTOKI installer/update.
	- `varargin`{: .text-blue-000} -- `(empty)`, `'demo'`, `'update'`.

```matlab
# USAGE
sorotoki();		% calls the installer
sorotoki('update');	% updates SOROTOKI to newest version
sorotoki('demo');	% provides a list of demos
```
</div>

<div class="code-example" markdown="1">
## `clr()`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Function`{: .text-red-000} --  Performs `clc`, `clear all`, and  `close all` 
</div>

<div class="code-example" markdown="1">
## `cdsoro()`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Function`{: .text-red-000} --  Sets current directory to installation folder of SOROTOKI.
</div>

<div class="code-example" markdown="1">
## `add2path()`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Function`{: .text-red-000} --  Adds current directory to MATLAB's search path.
</div>

<div class="code-example" markdown="1">
## `unload_sorotoki`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Function`{: .text-red-000} --  Removes the entire SOROTOKI toolkit from MATLAB's search path. **NOTICE:** This does not uninstall SOROTOKI, it prevents MATLAB from finding all functions tied to the toolkit. If MATLAB restarts, the `startup.m` will load SOROTOKI normally.
</div>


---

## Signed Distance Functions  -- `Sdf.m`{: .text-purple-000}
<div class="code-example" markdown="1">
## `sdf = Sdf(fnc)`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Class::Sdf`{: .text-red-000} --  Creates a Signed Distance Function Class based on `fnc = @(x) ....`
	- `fnc`{: .text-blue-000} -- `Function::f = @(x) ...`{: .text-red-000} such that `d = f([Nx2 Matrix])` or `d = f([Nx3 Matrix])` gives the output `d = [Nx1 Column]` of signed distances (negative implies inside the domain). The simplest example is `sdf = @(x) sqrt((x(:,1)).^2 + (x(:,2)).^2) - 1.0` which results in a unit-circle about the origin (0,0).

- `Public variables`{: .text-red-000} 
	- `sdf`{: .text-blue-000} -- `Function::sdf = @(x) ...`{: .text-red-000},
	- `BdBox`{: .text-blue-000} -- `[1x4 Row]`, `[1x6 Row]`,
	- `cmap`{: .text-blue-000} -- `viridis` (default), or `[Nx3 ColorMatrix]`.

```matlab
# USAGE
fnc = @(x) sqrt((x(:,1)).^2 + (x(:,2)).^2) - 1.0;
sdf = Sdf(fnc,'BdBox',[-1,1,-1,1]);
```
</div>

<div class="code-example" markdown="1">
## `Sdf = Sdf1 + Sdf2 + ... + Sdfn`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Class operator`{: .text-red-000} --  Unions two or more Sdf classes.
	- `Sdf1,Sdf2,...`{: .text-blue-000} -- `Class::Sdf`{: .text-red-000}
	- `Sdf`{: .text-blue-000} -- `Class::Sdf`{: .text-red-000}
</div>

<div class="code-example" markdown="1">
## `Sdf = Sdf1 - Sdf2 - ... - Sdfn`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Class operator`{: .text-red-000} --  Difference between two or more Sdf classes. `Sdf1` is the base function on which the operation is performed.
	- `Sdf1,Sdf2,...`{: .text-blue-000} -- `Class::Sdf`{: .text-red-000}
	- `Sdf`{: .text-blue-000} -- `Class::Sdf`{: .text-red-000}
</div>

<div class="code-example" markdown="1">
## `Sdf = Sdf1/Sdf2`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Class operator`{: .text-red-000} --  Intersection between two Sdf classes. `Sdf1` is the base function on which the operation is performed.
	- `Sdf1,Sdf2`{: .text-blue-000} -- `Class::Sdf`{: .text-red-000}
	- `Output`{: .text-blue-000} -- `Class::Sdf`{: .text-red-000}
</div>

<div class="code-example" markdown="1">
## `Sdf.show()`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Public function`{: .text-red-000} --  Creates `figure(101)` or uses existing `figure(101)` to show the Signed Distance Field within the domain `Sdf.BdBox`. The colormap is `viridis` by default.
</div>


---

## 2D/3D mesh generation -- `Mesh.m`{: .text-purple-000}

bla

---

## Finite element method -- `Fem.m`{: .text-purple-000}

bla

---


## Dynamic Modeling -- `Model.m`{: .text-purple-000}

bla

---


## Graphical models -- `Gmodel.m`{: .text-purple-000}

bla

---


## IK-rigging -- `Rig.m`{: .text-purple-000}

bla

---

## Plotting tools
<div class="code-example" markdown="1">
## `background(color)`{: .fs-6 .text-purple-000 .text-alpha}	
{: .no_toc }
- `Function`{: .text-red-000} --  Sets figures background color
	- `color`{: .text-blue-000} -- `'w'`, `'b'`, `gitepage`, `metropolis` 
</div>


---

## Linear- and hyper-elastic material presets

## Graphical material presets

bla

---