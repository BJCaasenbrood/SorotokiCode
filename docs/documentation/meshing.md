---
layout: default
title: Mesh generation
parent: Documentation
nav_order: 1
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 

<!-- <div align="center"> <img src="./src/mesh.png" width="650"> </div> -->

# Mesh generation
SOROTOKI offers mesh generation for triangular, quadrilateral, and polygonal elements. The restricted material domains are defined by so-called '*signed distance functions*' (SDF). The toolkit provides a set of these geometric shape SDFs (e.g., circles, rectangles, lines) and boolean operators (e.g., union, difference, and intersect). Together they allow for a wide range of spatial discretization in $$\mathbb{R}^2$$ and $$\mathbb{R}^3$$. 

# Mathematical background
## Signed distance functions
A signed distance functions (SDF) passes a spatial coordinate and returns the shortest distance to the boundaries of a metric domain. Mathematically speaking, the signed distance function $$d_\Omega: \mathbb{R}^n \mapsto \mathbb{R}$$ associated with the subset $$\Omega$$ of Euclidean space $$\mathbb{R}^n$$ is defined as

$$ d_\Omega(x) := s_\Omega(x) \min_{y \in \partial \Omega} \lVert x - y \rVert,$$ 

where $$s_\Omega$$ represents a discontinuous sign function

$$ s_\Omega(x) = 
\begin{cases}
-1, & x \in\Omega, \\
+1, & x \in \mathbb{R}^n\setminus \Omega,
\end{cases}
$$

and $$\partial \Omega$$ is the boundary of the material domain $$\Omega$$. The sign of the distance function determines if the coordinate is located inside or outside the bounded domain. Therefore, evaluation of the SDF function is not only numerically efficient, it also allows for an implicit representation of the spatial domain, which can be easily paired with various mathematical operations, like addition, subtraction, and differentiation. The toolkit comes with some preset SDF functions which can be used with off-the-shelf available mathematical operators.

### Preset SDF-functions
```rust
% two-dimensional
d = dCircle(P,xc,xy,r);
d = dLine(P,x1,x2,y1,y2);
d = dRectangle(P,x1,x2,y1,y2);

% three-dimensional
d = dCube(P,x1,x2,y1,y2,z1,z2);
d = dSphere(P,xc,yc,zc,r);
d = dCuboid(P,a,b,c);
```


### Mathematical operators

| Operation   | SOROTOKI     | Math  |
|:-------------|:------------------|:------|
| Build    | `OMEGA = dCircle(P,a,b,R)` | $$\Omega = \left\{r\in \mathbb{R}^2 \; : \; (r_x-a)^2 + (r_y - b)^2 \le R^2 \right\} $$ |
| Union    | `D = OMEGA1 + OMEGA2` | $$\mathcal{D} = \Omega_1 \cup \Omega_2$$  |
| Difference    | `D = OMEGA1 - OMEGA2` | $$\mathcal{D} = \Omega_1  \backslash \Omega_2$$  |
| Intersection    | `D = OMEGA1\OMEGA2` | $$\mathcal{D} = \Omega_1 \cap \Omega_2$$  |
| Evaluation    | `y = D.eval(x)` | $$y = d_\Omega(x)$$ with $$x\in \Omega$$ |


### Example

```rust
%% set signed distance function
msh = Mesh(@(x) SDF(x,0),'BdBox',[-1 3 -1 3]);
subplot(2,2,1); msh.showSDF();

msh = Mesh(@(x) SDF(x,1),'BdBox',[-2 2 -1 3]);
subplot(2,2,2); msh.showSDF();

msh = Mesh(@(x) SDF(x,2),'BdBox',[-1.5 2.5 -1 3]);
subplot(2,2,3); msh.showSDF();

msh = Mesh(@(x) SDF(x,3),'BdBox',[-1 3 -1 3]);
subplot(2,2,4); msh.showSDF();

function Dist = SDF(x,request)
    R = dRectangle(x,0,2,0,2);
    C = dCircle(x,0,1,1);

    switch(request)
        case(0), Dist = R;
        case(1), Dist = C;
        case(2), Dist = dUnion(R,C);
        case(3), Dist = dDiff(R,C);
    end
end
```