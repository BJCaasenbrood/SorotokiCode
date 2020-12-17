---
layout: default
title: Home
nav_order: 1
permalink: /
image_sliders:
  - example_slider
---


<!--
  https://fontawesome.com/how-to-use/on-the-web/styling/sizing-icons
-->

{% include slider_styles.html %}

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 
<div align="center"> <img src="./docs/documentation/img/softrobot_.png" width="650"> </div> <br/>

# **SOROTOKI** - An open-source soft robotics toolkit for MATLAB
{: .fs-6 .text-purple-000}

*SOROTOKI* is an open-source MATLAB toolkit for soft robotics that includes an array tools for design, modeling, and control. Due to its scientific diversity, it can be challenging for researchers to quickly learn the different disciplines required for soft robotics. With the aim to lower the multi-disciplinary barrier, Sorotoki aims to incorporate multiple layers of scientific disciplines in one toolkit. Examples include: continuum mechanics, dynamic system- and control theory, topology optimization, computer graphics, and much more to come! 

## List of applications:
{: .text-purple-000}
- [x] Modeling of 2D/3D geometries with Signed Distance Functions (SDFs),
- [x] Performing Finite Element Analysis (FEA) using hyper-elastic materials,
- [x] Topology optimization-based design of (pressure-driven) soft robots,
- [x] Dynamic modeling of soft robots through Cosserat-beam theory,
- [x] Fast and responsive rendering with automatically-updated textures.

#### To do:
{: .text-purple-000}
- [ ] Real-time control of physical soft robots via BalloonDog interface,
- [ ] Including class for magnetic fields/inductive sensing,
- [ ] A pdf-version of the documentation.

{% include slider.html selector="example_slider" %}

### Download the latest stable version
{: .text-purple-000}
[Stable V3.03 (.zip)](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 .mr-2} [Stable V3.03 (.tar)](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master){: .btn .btn-purple .fs-5 .mb-4 .mb-md-0 .mr-2} [View on Github](https://github.com/BJCaasenbrood/SorotokiCode){: .btn .fs-5 .mb-4 .mb-md-0}  

or via git clone
{: .fs-5 .text-purple-000}
```fortran
 git clone https://github.com/BJCaasenbrood/SorotokiCode.git
```

### Citation
{: .text-purple-000}
If you are using Sorotoki in your (academic) work, please consider citing the toolkit:
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


<!-- ## Quick installation
The toolkit is easy to install. Download the latest stable version above, and unpack the compressed folder at any desired directory. Alternatively, you can directly clone the repository with the following terminal command:
```rust
git clone https://github.com/BJCaasenbrood/SorotokiCode.git
```

To install the toolkit, simply run the command [`sorotoki()`]() in the MATLAB command window. That's it, the SOROTOKI toolkit is now ready-to-use.

It is worth mentioning this command is also used to update the toolkit. It is recommended to run `sorotoki.m`{: .text-purple-000} to check for updates occasionally.  -->

<!-- ## Citation and references
If you are using Sorotoki in your (academic) work, please consider citing the toolkit:
```bibtex
@misc{Caasenbrood2018,
  author = {Caasenbrood, Brandon},
  title = {Sorotoki - A Soft Robotics Toolkit for MATLAB},
  year = {2018},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/BJCaasenbrood/SorotokiCode}},
}
```

[**[1]** ](https://ieeexplore.ieee.org/abstract/document/9116010/metrics#metrics) B. Caasenbrood, A. Pogromsky and H. Nijmeijer, A Computational Design Framework for Pressure-driven Soft Robots through Nonlinear Topology Optimization, 2020 3rd IEEE Inter. Conf. on Soft Robotics (RoboSoft), pp. 633-638, 2020.

[2] B. Caasenbrood, A. Pogromsky, and H. Nijmeijer, Dynamic modeling of hyper-elastic soft robots using spatial curves, IFAC World Congress, 2020.

[**[3]**](https://link.springer.com/article/10.1007/s00158-011-0706-z) C. Talischi, G. H. Paulino, A. Pereira, and I. F. M. Menezes, PolyMesher: A general-purpose mesh generator for polygonal elements written in Matlab, Struct. Multidiscip. Optim., vol. 45, no. 3, pp. 309–328, 2012.

[**[4]**](https://www.springer.com/gp/book/9781441917454) N. Kim, Introduction Analysis Finite Element to Nonlinear. Springer, 2018.

[**[5]**](https://www.springer.com/gp/book/9783540429920) M. Bendsøe, M. Philip, O. Sigmund, Topology Optimization. Theory, Methods and Applications. Springer, 2003.

 -->
