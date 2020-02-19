<div align="center"> <img src="./src/fem.png" width="650"> </div>

# Nonlinear Finite Elements

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)

### List of material models:
```matlab
% different material models
fem.Material = LinearMaterial('E',-,'Nu',-);           % linear isotropic
fem.Material = NeoHookeanMaterial('E',-,'Nu',-);       % neo-hookean model
fem.Material = MooneyMaterial('C10',-,'C01',-,'K',-);  % mooney-rivlin model
fem.Material = YeohMaterial('C1',-,'C2',-,'C3',-)      % yeoh model
                            
% different preset materials
fem.Material = Ecoflex0030();
fem.Material = Dragonskin20A();   	 
fem.Material = Elastosil();   	 
```

### Example: Clamped beam 
<div align="center"> <img src="./src/fem_beam.png" width="350"> </div>

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
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[0,-1e-3]);

%% assign material
fem.Material = Ecoflex0030;

%% solving
fem.solve();
```

### Example: Tensile bone
```matlab
%% generate mesh from sdf
sdf = @(x) TensileBone(x,8,2,3,1,0.75);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,10,0,10],'NElem',500);
msh = msh.generateMesh;

%% generate fem model from mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/10,'PrescribedDisplacement',true);

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,4]);

%% assign material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

function D = TensileBone(P,H,W,T,D,R)
dD = 0.5*(W-D);
dT = 0.5*(H-T);

R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,0,dD,dT,dT+T);
R3 = dRectangle(P,W-dD,W,dT,dT+T);
C1 = dCircle(P,dD-R,dT,R);
C2 = dCircle(P,dD-R,dT+T,R);
C3 = dCircle(P,W-dD+R,dT,R);
C4 = dCircle(P,W-dD+R,dT+T,R);
D0 = dDiff(dDiff(dDiff(R1,R2),C1),C2);
D = dDiff(dDiff(dDiff(D0,R3),C3),C4);
end
```
