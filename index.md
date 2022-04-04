---
layout: default
title: Home
nav_order: 1
permalink: /
---
<!--
  https://fontawesome.com/how-to-use/on-the-web/styling/sizing-icons
-->

<!-- <script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> -->

<div align="center"> <img src="./docs/documentation/img/softrobot_.png" width="650"> </div> <br/>

# **SOROTOKI** - An open-source soft robotics toolkit for MATLAB
{: .fs-6 .text-purple-000}

*SOROTOKI* is an open-source MATLAB toolkit for soft robotics that includes an array of tools for design, modeling, and control. Due to its scientific diversity, it can be challenging for new researchers to quickly familiarize themselves with multiple scientific areas. With the aim to lower this threshold, Sorotoki aims to incorporate multiple layers of soft robotics research into one compact toolkit. Examples include: continuum mechanics, dynamic systems and control theory, topology optimization, computer graphics, and much more to come! The combination provides a highly flexible programming environment and will hopefully aid the development of novel soft robotic research.

#### Download the latest stable version:
{: .text-purple-000}

```bash
git clone --depth 1  https://github.com/BJCaasenbrood/SorotokiCode.git
```

{: .text-purple-000}
[Stable V2.24.01  (.zip)](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 .mr-2} [Stable V2.24.01 (.tar)](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 .mr-2} [View on Github](https://github.com/BJCaasenbrood/SorotokiCode){: .btn .fs-5 .mb-4 .mb-md-0}  

## What's new?
{: .text-purple-000}
- **Apr 03, 2022**: We added environmental contact in SOROTOKI! You can now use Signed Distance Functions `Sdf.m`{: .text-purple-000} with `Fem.m`{: .text-purple-000}. Use `fem.AddConstraint('Contact',sdf,Move)`{: .text-purple-000} to add the contact. The vector `Move = [dX, dY]`{: .text-purple-000} or `Move = [dX, dY, dZ]`{: .text-purple-000} can be used to offset the contact. The contact dynamics can be solved Quasi-static or Dynamic.

<details>
<summary>Show image</summary>
<div align="center"> <img src="./docs/documentation/img/fem_bouncing_ball.gif" width="500"> </div>
</details>

- **Mar 10, 2022**: We hugely improved the dynamic simulations performance using a combination of Matrix-Differential Equations (MDEs) and Matlab's mex compiler. To enable the compiled dynamic solver, set `mdl.MexSolver = true`{: .text-purple-000} (default option). As an example, see `mdl_swing_motion`{: .text-purple-000}. It runs a 16-dim Kirchoff beam at 60 Hz easily! (5s simulation = 1s real-time).
Specs: Lenovo Yoga, AMD Ryzen 7 5800H, 32GB RAM.

<details>
<summary>Show video</summary>
<div align="center"> <img src="./docs/documentation/img/mdl_fast_simu.mp4" width="500"> </div>
</details>

<a href="https://github.com/BJCaasenbrood/SorotokiCode/blob/master/scripts/mdl/mdl_swing_motion.m"> Code available here</a>

- **Nov 11, 2021**: We recently added the functionality to perform real-time control of
soft robotic systems using a desktop-sized development platform, see [Github repo](https://github.com/chukhanhhoang/SoftRoboticSetupFesto). Using `Model.m`{: .text-purple-000}, we developed a model-based controller which we see in action below!

<details>
<summary>Show image</summary>
<div align="center"> <img src="./docs/documentation/img/bdog_closed_loop_control.gif" width="500"> </div>
</details>

{: .text-purple-000}
**REMARK**: The controller shown in the video is not a 'elementary' PID controller, we actively use gravity compensation, and hyper-elastic and visco-elastic material compensation! To envelop these entities in a dynamic model, we used the SOROTOKI toolkit.
{: .fs-2}


## Applications highlights
{: .text-purple-000}
- Implicit modeling with Signed Distance Functions (SDFs),
- Finite element method (FEM) using hyper-elastic materials,
- Topology optimization of (pressure-driven) soft robots,
- Dynamical modeling through differential geometric theory,
- (NEW!) Real-time control of soft robots via Raspi-interface,
- Fast graphics rendering with responsive textures.

{: .text-purple-000}
**REMARK**: All images below are produced using only MATLAB and the SOROTOKI toolkit, no additional software was used!
{: .fs-2}

### Signed Distance Functions and Meshing -- `Sdf.m`{: .text-purple-000}, `Mesh.m`{: .text-purple-000}

#### SDF: Implicit modeling of 2D-primes: union and subtraction
{: .text-purple-000}
<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/sdf_primes.png" width="500"> </div>

<a href="https://github.com/BJCaasenbrood/SorotokiCode/blob/master/scripts/mesh/mesh_sdf.m">Code available here</a>
</details>

#### SDF: Implicit modeling of 3D-primes: subtraction and intersection
{: .text-purple-000}
<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/sdf_3Dprimes.png" width="250"> </div>

<a href="https://github.com/BJCaasenbrood/SorotokiCode/blob/master/scripts/gmdl/SDF/preview_sdf.m">Code available here</a>
</details>

#### MESH: Signed distance function (SDF) to mesh
{: .text-purple-000}
<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/msh_lens.png" width="500"> </div>
</details>

#### MESH: Polygonal meshing of circular SDF
{: .text-purple-000}
<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/msh_circle_gen.gif" width="500"> </div>
</details>

### Finite Element Method  -- `Fem.m`{: .text-purple-000}

#### FEM: Uni-axial tensile test
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/fem_tensile.gif" width="500"> </div>

<a href="https://github.com/BJCaasenbrood/SorotokiCode/blob/master/scripts/fem/2D/static/fem_tbone.m">Code available here</a>
</details>

#### FEM: Topology optimization of PneuNet actuator
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/opt_pneunet.gif" width="500"> </div>

<a href="https://github.com/BJCaasenbrood/SorotokiCode/blob/master/scripts/fem/2D/static/fem_tbone.m">Code available here</a>
</details>

#### FEM: Deforming PneuNet actuator -- Ecoflex 0030
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/fem_pneunet.gif" width="500"> </div>
</details>

#### FEM: 3D buckling of hyper-elastic beam
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/fem_3D_buckle.gif" width="200"> </div>
</details>

#### FEM: Hyper-elastic bouncing ball
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/fem_bouncing_ball.gif" width="400"> </div>
</details>

#### FEM: PneuNet with dynamic contact
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/fem_pneunet_dyn_contact.gif" width="300"> </div>
</details>


### Dynamic Model  -- `Model.m`{: .text-purple-000}

#### MODEL: Simulation of two-link soft robot
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/mdl_softarm.gif" width="500"> </div>
</details>

#### MODEL: Energy-based control of planar soft robot
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/mdl_twolink_control.gif" width="500"> </div>
</details>

#### MODEL: Basis reconstruction from FEM-data -- Reduced-order modeling
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/mdl_femconstruct.png" width="900"> </div>
</details>

### Real-time Control -- `Control.m`{: .text-purple-000}

#### CONTROL: Closed-loop control of PneuNet actuator
{: .text-purple-000}
<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/bdog_closed_loop_control.gif" width="400"> </div>
</details>

#### CONTROL: Energy-based control two-link system
{: .text-purple-000}
<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/bdog_2link.gif" width="300"> </div>
</details>

### Graphics Model -- `Gmodel.m`{: .text-purple-000}

#### GMODEL: Responsive rendering of the Stanford bunny
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/examples/graphics/renderstl/img/rotateBunny.gif" width="500"> </div>
</details>

#### GMODEL: Responsive lighting
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/david.gif" width="300"> </div>
</details>

#### GMODEL: Rendering ambient occlusion (AO)
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/bunny_AO.png" width="550"> </div>
</details>

#### GMODEL: Rendering sub-surface scattering (SSS)
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/bunny_SSS.png" width="550"> </div>
</details>

## Citation
{: .text-purple-000}

If you are planning on using Sorotoki in your (academic) work, please consider citing the toolkit  

```bibtex
@misc{Caasenbrood2018,
  author = {Caasenbrood, Brandon},
  title = {Sorotoki - A Soft Robotics Toolkit for MATLAB},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/BJCaasenbrood/SorotokiCode}},
}
```
