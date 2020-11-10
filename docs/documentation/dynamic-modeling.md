---
layout: default
title: Dynamical modeling
parent: Documentation
nav_order: 4
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 

# Dynamic modeling 
[**Sorotoki**](https://bjcaasenbrood.github.io/SorotokiCode/) **Still under developement!** We have some interesting simulations though.



$$g(\sigma,t) := \begin{pmatrix} \Phi(\sigma,t) & r(\sigma,t) \\ 0_3^\top & 1 \end{pmatrix} \in \mathbb{SE}(3),
$$

where $$\Phi \in \mathbb{SO}(3)$$ is an orientation matrix, and $$r\in \mathbb{R}^3$$ is a position vector. Following the differential geometry of the Lie groups, the partial derivates of $$g \in \mathbb{SE}(3)$$ w.r.t. time and space can be expressed by two vector fields belonging to the Lie algebra of $$\mathbb{SE}(3)$$, which are defined by

$$\frac{\partial g}{\partial t} = g \hat{\eta} \quad \Longrightarrow \quad \eta(\sigma,t) := \begin{pmatrix} \Omega(\sigma,t) \\  V(\sigma,t) \end{pmatrix},
$$

and

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
