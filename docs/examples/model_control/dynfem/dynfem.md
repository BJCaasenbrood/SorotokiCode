---
layout: default
title: Example 3 - Reduced-order Models using FEM-driven Galerkin Projections.
parent: Modeling and Control 
grand_parent: Examples
nav_order: 3
permalink: /docs/examples/model_control/dynfem
---

#  Reduced-order dynamics using FEM-driven Galerkin projection (snapshot method).
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

#### Difficulty: `hard`{: .fs-3 .text-red-200} 
{: .no_toc }
 - Required classes: `Mesh.m`{: .text-purple-000}, `Fem.m`{: .text-purple-000}
 - Code length: `~25 lines`{: .text-purple-000} (without comments)

---

### Introduction
In this illustrative example, we will perform a simple pull test using a hyper-elastic material -- Ecoflex-0030 from SmoothOn. Assuming a two-dimensional problem, we consider a 20x20 specimen and perform an uni-axial elongation of $\lambda_{1} = 500\%$ (with $\lambda_2 = \lambda_3$). We model this using a single quadrilateral finite-element subjected to plane-stress conditions.

<div align="center"> <img src="./img/result.png" width="500"> </div>
<div align="center"> Stock image of uni-axial test (left, see [1]). Produced result from SOROTOKI (right). </div>

[**[1]**](https://www.rubbernews.com/blogs/new-products-james-heal-makes-move-rubber-testing-markets) **Titan10** a pull-test instrument for rubber and elastomers.
{: .fs-3} 

