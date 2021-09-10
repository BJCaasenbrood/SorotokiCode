---
layout: default
title: About SOROTOKI
parent: Documentation
nav_order: 1
---

# About
SOROTOKI is an open-source MATLAB toolkit for soft robotics that includes an array of tools for design, modeling, and control. Due to its scientific diversity, it can be challenging for researchers to quickly familiarize themselves with multiple scientific disciplines. With the aim to lower the threshold, SOROTOKI aims to incorporate multiple layers of soft robotics research into one toolkit. Examples include: continuum mechanics, dynamic systems and control theory, topology optimization, computer graphics, and much more to come! The combination provides a highly flexible modeling environment and hopefully aids the development of novel soft robotic research.

# Object-oriented architecture
SOROTOKI is an object-oriented toolkit that aims for a minimalistic code style. To do so, the toolkit is equipped with several Classes with user-input functionalities, which allow for cross-interaction between other Classes. The main classes are **Sdf**, **Mesh**, **Fem**, **Model**, **Gmodel**, and **Bdog**.  A flow-diagram of the toolkit on how these Classes work and interact is given schematically below:

<div align="center"> <img src="./img/diagram.svg" height="1520"> </div>

### The main object-oriented classes

#### Signed-distance function -- `Sdf.m`

#### Mesh generation -- `Mesh.m`
