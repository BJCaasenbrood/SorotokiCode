<div align="center"> <img src="./src/fem.png" width="650"> </div>

# Nonlinear Topology Optimization

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)

### Example: Beam 
```matlab
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,10,0,2);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,2],'NElem',500);
msh = msh.generateMesh;

%% generate fem model from mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/15,'ResidualNorm',1e-3);

%% add boundary conditions 
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[0,-2e-3]);

%% assign material
fem.Material = Ecoflex0030;

%% solving
fem.solve();
```

