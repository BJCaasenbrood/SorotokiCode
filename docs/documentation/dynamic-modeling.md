---
layout: default
title: Dynamical modeling
parent: Documentation
nav_order: 4
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 

# Dynamic modeling 
[**Sorotoki**](https://bjcaasenbrood.github.io/SorotokiCode/) **Still under developement!** We have some interesting simulations though.


# Cosserat beam theory 
Let $$\sigma \in \mathbb{X}$$ be a spatial variable that lies on a bounded domain $$\mathbb{X} \in [0,l]$$ with $$l \in \mathbb{R}_+ $$ the intrinsic length of the soft body, and let $$t \in \mathbb{R}$$ be a temporal variable. Given this coordinate system, we can represent any rigid-body frame that is rigidly attached to the deformable body of the soft robot by a space-time variant curve $$g : \mathbb{X} \times \mathbb{R} \mapsto \mathbb{SE}(3)$$, that is,

$$g(\sigma,t) := \begin{pmatrix} \Phi(\sigma,t) & r(\sigma,t) \\ 0_3^\top & 1 \end{pmatrix} \in \mathbb{SE}(3),
$$

where $$\Phi \in \mathbb{SO}(3)$$ is an orientation matrix, and $$r\in \mathbb{R}^3$$ is a position vector. Here, the $$\mathbb{SE}(3)$$ is a Lie group associated with rigid-body transformations on $$\mathbb{R}^3$$, whose smoothness allows differentiable operations through its respective Lie algebra. Following the differential geometry of the Lie groups, the partial derivates of $$g \in \mathbb{SE}(3)$$ w.r.t. time and space can be defined by two vector fields belonging to the Lie algebra of $$\mathbb{SE}(3)$$, which are defined by

$$\frac{\partial g}{\partial t} = g \hat{\eta} \quad \Longrightarrow \quad \eta(\sigma,t) := \begin{pmatrix} \Omega(\sigma,t) \\  V(\sigma,t) \end{pmatrix},
$$

$$\frac{\partial g}{\partial \sigma} = g \hat{\xi} \quad \Longrightarrow \quad \xi(\sigma,t) := \begin{pmatrix} K(\sigma,t) \\  U(\sigma,t) \end{pmatrix},
$$

respectively. 

The total (co)-energy function or Hamiltonian of the continuous dynamical system is then given by

$$ \mathcal{H} = \int_\mathbb{X} \eta^\top \mathcal{M} \eta \; d\sigma + \int_\mathbb{X} W(\xi) \; d\sigma$$

## Soft robot manipulator dynamic simulation
<div align="center"> <img src="./img/SoftArm.gif" width="500"> </div>

## Pneu-net spatial tracking (model-based controller)
<div align="center"> <img src="./img/Pneunet_tracking.gif" width="500"> </div>

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)
