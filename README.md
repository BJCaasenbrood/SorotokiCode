<div align="center"> <img src="https://raw.githubusercontent.com/BJCaasenbrood/SorotokiCode/gh-pages/img/softrobot.svg" width="850"> </div> <br/>

[![CircleCI](https://circleci.com/gh/BJCaasenbrood/SorotokiCode.svg?style=svg)](https://circleci.com/gh/BJCaasenbrood/SorotokiCode) [![MATLAB](https://github.com/BJCaasenbrood/SorotokiCode/actions/workflows/CI_Linux.yml/badge.svg)](https://github.com/BJCaasenbrood/SorotokiCode/actions/workflows/CI_Linux.yml) [![Static Badge](https://img.shields.io/badge/Sorotoki-documentation-blue)](https://bjcaasenbrood.github.io/SorotokiCode/) [![Discord](https://img.shields.io/discord/1049668116688412672?logo=discord)](https://discord.gg/u6CWsjPsyG) ![GitHub](https://img.shields.io/github/license/BJCaasenbrood/SorotokiCode) 

[**Sorotoki**](https://bjcaasenbrood.github.io/SorotokiPage/) is an open-source MATLAB toolkit for soft robotics that includes an array of tools for design, modeling, and control. Due to its scientific diversity, it can be challenging for new researchers to quickly familiarize themselves with multiple scientific areas. With the aim to lower this threshold, Sorotoki aims to incorporate multiple layers of soft robotics research into one compact toolkit. Examples include: continuum mechanics, dynamic systems and control theory, topology optimization, computer graphics, and much more to come! The combination provides a highly flexible programming environment and will hopefully aid the development of novel soft robotic research.

#### Installation via Matlab Package Installer (MPI)
The toolkit is easy to install using the Matlab Package Installer ([mpi](https://github.com/mobeets/mpm)). Please follow the [**Installation Guide**](https://bjcaasenbrood.github.io/SorotokiCode/install/).


## Applications highlights

- [x] [Sdf](https://github.com/BJCaasenbrood/SorotokiSdf) - Implicit modeling with Signed Distance Functions (SDFs),
- [x] [Fem](https://github.com/BJCaasenbrood/SorotokiFem) - Finite element method (FEM) using hyper-elastic materials,
- [x] [Topo](https://github.com/BJCaasenbrood/SorotokiFem) - Topology optimization of (pressure-driven) soft robots,
- [x] [Model](https://github.com/BJCaasenbrood/SorotokiModel) - Dynamical modeling through geometric theory,
- [x] [Control](https://github.com/BJCaasenbrood/SorotokiCode) - Real-time control of soft robots via Raspi-interface,
- [x] [GModel](https://github.com/BJCaasenbrood/MatlabGraphicalModel) - Fast graphics rendering with responsive textures.

## What can it do?

<img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/mckibben.gif" height="150">  <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/twist.gif" height="140"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/buckling.png" height="140"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/diamondbot.png" height="120"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/opt_bellow.gif" height="150"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/opt_pneunet_90.gif" height="160">

<img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/straingauge.gif" height="145"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/soft_finger.gif" height="140"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/soro_softcrawl.gif" height="140"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/soft_bounce.gif" height="140"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/soro_gripper.png" height="140">


<img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/beam.gif" height="160"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/soro_hand.gif" height="150"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/soft_control.gif" height="150"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/soro_control.gif" height="140"> 

<img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/bdog_bellow.gif" height="160"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/bdog_grip.gif" height="160"> <img src="https://github.com/BJCaasenbrood/SorotokiPage/blob/master/docs/documentation/img/bdog_srm.gif" height="160"> 

## Citation

If you are planning on using Sorotoki in your (academic) work, please consider citing the toolkit  

```bibtex
@article{Caasenbrood2024,
	author = {Caasenbrood, Brandon J. and Pogromsky, Alexander Y. and Nijmeijer, Henk},
	title = {{Sorotoki: A Matlab Toolkit for Design, Modeling, and Control of Soft Robots}},
	journal = {IEEE Access},
	volume = {12},
	pages = {17604--17638},
	year = {2024},
	publisher = {IEEE},
	doi = {10.1109/ACCESS.2024.3357351}
}
```
