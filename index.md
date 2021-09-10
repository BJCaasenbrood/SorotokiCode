---
layout: default
title: Home
nav_order: 1
permalink: /
---
<!--
  https://fontawesome.com/how-to-use/on-the-web/styling/sizing-icons
-->

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
<div align="center"> <img src="./docs/documentation/img/softrobot_.png" width="650"> </div> <br/>

# **SOROTOKI** - An open-source soft robotics toolkit for MATLAB
{: .fs-6 .text-purple-000}

*SOROTOKI* is an open-source MATLAB toolkit for soft robotics that includes an array of tools for design, modeling, and control. Due to its scientific diversity, it can be challenging for new researchers to quickly familiarize themselves with multiple scientific areas. With the aim to lower this threshold, Sorotoki aims to incorporate multiple layers of soft robotics research into one compact toolkit. Examples include: continuum mechanics, dynamic systems and control theory, topology optimization, computer graphics, and much more to come! The combination provides a highly flexible programming environment and will hopefully aid the development of novel soft robotic research.

#### Download the latest stable version (v1.03):
{: .text-purple-000}

```bash
git clone --depth 1  https://github.com/BJCaasenbrood/SorotokiCode.git
```

{: .text-purple-000}
[Stable V1.03 (.zip)](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 .mr-2} [Stable V1.03 (.tar)](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 .mr-2} [View on Github](https://github.com/BJCaasenbrood/SorotokiCode){: .btn .fs-5 .mb-4 .mb-md-0}  

## Applications highlights
{: .text-purple-000}
- [x] Implicit modeling with Signed Distance Functions (SDFs),
- [x] Finite element method (FEM) using hyper-elastic materials,
- [x] Topology optimization of (pressure-driven) soft robots,
- [x] Dynamical modeling through geometric theory,
- [x] (NEW!) Real-time control of soft robots via Raspi-interface,
- [x] Fast graphics rendering with responsive textures.

{: .text-purple-000}

**REMARK**: All images are produced using only MATLAB and SOROTOKI, no additional software was used!
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
</details>

#### FEM: Topology optimization of PneuNet actuator
{: .text-purple-000}

<details>
<summary>Show image</summary>
<div align="left"> <img src="./docs/documentation/img/opt_pneunet.gif" width="500"> </div>
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
