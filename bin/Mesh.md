<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 
<div align="center"> <img src="./src/mesh.png" width="650"> </div>

# Mesh Generation
[**Sorotoki**](https://bjcaasenbrood.github.io/SorotokiCode/) offers mesh generation for triangular, quadrilateral, and polygonal elements. The restricted material domains for the meshes are defined by so-called *signed distance functions* (SDF). The toolkit provides a set of geometeric shape (e.g., circles, rectangles, lines) and boolean operators, e.g., union, difference, and intersect. Together these operation allow for a wide range of meshing domain. 

# Signed distance functions
A signed distance functions (SDF) passes spatial coordinates and returns a signed value representing the shortest distance to the boundary of a spatial domain $\Omega$. Mathematically, the signed distance function $$f: \mathbb{R^n} \mapsto \mathbb{R}$$ of a subset $$\omega \subset \R^n$$ is defined by

$$
    f(x) := 
\begin{cases}
    d(x,\partial \Omega )   & \text{if } x \in \Omega\\
    -d(x,\partial \Omega ) & \text{if } x \in \Omega/\R^n
\end{cases}
$$

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)
