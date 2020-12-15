---
layout: default
title: Dynamical modeling
parent: Documentation
nav_order: 4
---


<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 

# Dynamic modeling 
**Still under developement!** We have some interesting simulations though.


**Disclaimer:**{: .text-purple-000} To modify the embedded dynamics code of the SOROTOKI toolkit, you will required a c++ compiler (e.g., gcc/g++/clang) and the EIGEN linear algebra library for c++. It is not necessary to use CMake, however, it is highly recommended.

# Mathematical background
## Cosserat beam theory 
Cosserat beam theory is based on the assumption that slender structures can be modeled as a one-dimensional PDE model subjected to various three-dimensional forces. Unlike classical Finite-Element Methods (FEM), Cosserat beam theory only requires spatial discretization in one dimension -- along the abscissa of the beam. Therefore, it requires less computational effort than FEM and, more importantly, it lends itself for Lie group theory that is popularized recently in modern robotics. 

To elaborate briefly on the mathematics behind the Cosserat beam theory, consider the following. Let $$\sigma \in \mathbb{X}$$ be a spatial variable that lies on a bounded domain $$\mathbb{X} \in [0,l]$$ with $$l \in \mathbb{R}+ $$ the intrinsic length of the soft body, and let $$t \in \mathbb{R}$$ be the time. As such, we can represent any orientation frame that is rigidly attached to the continuously deformable body by a one-dimensional curve $$g : \mathbb{X} \times \mathbb{R} \mapsto \mathbb{SE}(3)$$, that is,

$$
g(\sigma,t) := \begin{pmatrix} \Phi(\sigma,t) & r(\sigma,t) \\ 0_3^\top & 1 \end{pmatrix} \in \mathbb{SE}(3),
$$

where $$\Phi \in \mathbb{SO}(3)$$ is an orientation matrix, and $$r\in \mathbb{R}^3$$ is a position vector. Here, the $$\mathbb{SE}(3)$$ is a Lie group associated with rigid-body transformations on $$\mathbb{R}^3$$, whose smoothness allows differentiable operations through its respective Lie algebra. 

## Relating deformation and velocity
Following the differential geometry of the Lie groups, the partial derivates of $$g \in \mathbb{SE}(3)$$ w.r.t. space and time are described by two vector fields belonging to its Lie algebra $$\mathbb{se}(3)$$. Consequently, the deformation twist and the velocity twist are given as follows

$$\tfrac{\partial g}{\partial \sigma} = g \hat{\xi} \quad \Longrightarrow \quad \xi(\sigma,t) := \begin{pmatrix} K(\sigma,t) \\  U(\sigma,t) \end{pmatrix};
$$

$$\tfrac{\partial g}{\partial t} = g \hat{\eta} \quad \Longrightarrow \quad \eta(\sigma,t) := \begin{pmatrix} \Omega(\sigma,t) \\  V(\sigma,t) \end{pmatrix};
$$

Here $$K$$ and $$U$$ are the twist/curvature and elongation/shear deformations, and $$\Omega$$ and $$V$$ are the angular and linear velocities, respectively. These two vector field are the main premise behind the equation of motion of the Cosserat beam. Using the commutative property of the Lie algebras, we can relate the velocity-twist to the deformation-twist through the following PDE:

$$ \tfrac{\partial \eta}{\partial \sigma}= -\text{ad}_{\xi}(\eta) + \tfrac{\partial \xi}{\partial t}$$

## Port-Hamiltonian systems
The total (co)-energy function or Hamiltonian of the continuous dynamical system is then given by

$$ \mathcal{H} = \int_\mathbb{X} \eta^\top \mathcal{M} \eta \; d\sigma + \int_\mathbb{X} \Psi(\xi) \; d\sigma$$

# Numerical examples
**Still under developement!** We have some interesting simulations though.
## Soft robot manipulator dynamic simulation
<div align="center"> <img src="./img/SoftArm.gif" width="500"> </div>

## Pneu-net spatial tracking (model-based controller)
<div align="center"> <img src="./img/Pneunet_tracking.gif" width="500"> </div>

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)
