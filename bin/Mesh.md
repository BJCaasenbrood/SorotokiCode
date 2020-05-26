<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 
<div align="center"> <img src="./src/mesh.png" width="650"> </div>

# Mesh Generation
[**Sorotoki**](https://bjcaasenbrood.github.io/SorotokiCode/) offers mesh generation for triangular, quadrilateral, and polygonal elements. The restricted material domains for the meshes are defined by so-called *signed distance functions* (SDF). The toolkit provides a set of geometeric shape (e.g., circles, rectangles, lines) and boolean operators, e.g., union, difference, and intersect. Together these operation allow for a wide range of meshing domain. 

# Signed distance functions
A signed distance functions (SDF) passes a spatial coordinate and returns the shortest distance to the boundary of a domain $$\Omega$$. In mathematical terms, the signed distance function $$d_\Omega: \mathbb{R}^n \mapsto \mathbb{R}$$ assosciated wit the subset $$\Omega$$ of Euclidean space $$\mathbb{R}^n$$ is defined by

$$ d(x) := s_\Omega(x) \min_{y \in \partial \Omega} \lvert x - y \rvert$$ 
$$ s_\Omega
\begin{cases}
-1, & x \in\Omega \\
+1, & x \in \R^n\setminus \Omega
\end{cases}
$$

where $$s_\Omega(x)$$ representing a sign function, and $$\partial \Omega$$ the boundary of the domain $$\Omega$$. The sign of the distance function determines if the coordinate is inside or outside the bounded domain. 

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)
