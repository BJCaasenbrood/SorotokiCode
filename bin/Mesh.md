<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 
<div align="center"> <img src="./src/mesh.png" width="650"> </div>

# Mesh Generation
[**Sorotoki**](https://bjcaasenbrood.github.io/SorotokiCode/) offers mesh generation for triangular, quadrilateral, and polygonal elements. The restricted material domains for the meshes are defined by so-called *signed distance functions* (SDF). The toolkit provides a set of geometeric shape (e.g., circles, rectangles, lines) and boolean operators, e.g., union, difference, and intersect. Together these operation allow for a wide range of meshing domain. 

# Signed distance functions
A signed distance functions (SDF) passes a spatial coordinate and returns the shortest distance to the boundary of a metric domain. Mathematically, the signed distance function $$d_\Omega: \mathbb{R}^n \mapsto \mathbb{R}$$ assosciated with the subset $$\Omega$$ of Euclidean space $$\mathbb{R}^n$$ is defined by

$$ d_\Omega(x) := s_\Omega(x) \min_{y \in \partial \Omega} \lVert x - y \rVert$$ 

$$ s_\Omega = 
\begin{cases}
-1, & x \in\Omega \\
+1, & x \in \mathbb{R}^n\setminus \Omega
\end{cases}
$$

where $$s_\Omega(x)$$ representing a sign function, and $$\partial \Omega$$ the boundary of the domain $$\Omega$$. The sign of the distance function determines if the coordinate is inside or outside the bounded domain. 

### Preset SDF
```matlab
% two-dimensional
d = dCircle(P,xc,xy,r);
d = dLine(P,x1,x2,y1,y2);
d = dRectangle(P,x1,x2,y1,y2);

% three-dimensional
d = dCube(P,x1,x2,y1,y2,z1,z2);
d = dSphere(P,xc,yc,zc,r);
d = dCuboid(P,a,b,c);

```

### Example

```matlab
%% set signed distance function
msh = Mesh(@(x) SDF(x,0),'BdBox',[-1 3 -1 3]);
subplot(2,2,1); msh.showSDF;

msh = Mesh(@(x) SDF(x,1),'BdBox',[-2 2 -1 3]);
subplot(2,2,2); msh.showSDF;

msh = Mesh(@(x) SDF(x,2),'BdBox',[-1.5 2.5 -1 3]);
subplot(2,2,3); msh.showSDF;

msh = Mesh(@(x) SDF(x,3),'BdBox',[-1 3 -1 3]);
subplot(2,2,4); msh.showSDF;

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

[**Top of page**](https://bjcaasenbrood.github.io/SorotokiCode/bin/Mesh.html)

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)
