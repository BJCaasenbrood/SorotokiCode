<div align="center"> <img src="./src/fem.png" width="650"> </div>

# Nonlinear Topology Optimization

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)

### Example: Pneu-net soft robot

```matlab
%% generate mesh from sdf
sdf = @(x) PneuNet(x,20,40,1,20);

msh = Mesh(sdf'BdBox',[0,20,0,40],'NElem',1000);
msh = msh.generateMesh;

%% generate fem model from mesh
fem = Fem(msh,'VolumeInfill',0.3,...
              'Penal',4,...
              'FilterRadius',4,...
              'OptimizationProblem','Compliant');
              
%% set spatial settings
fem = fem.set('Periodic',[1/2, 0],'Repeat',ones(1,7));

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,1]);

id = fem.FindElements('Location',[W/2,H/2],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[W/2,H/2],1);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();

function Dist = PneuNet(P,W,H,E,T)
R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,-W/2,E,T,H+H/2);
R3 = dRectangle(P,W-E,W+W/2,T,H+H/2);
C1 = dCircle(P,0,T + 0.5,1);
C2 = dCircle(P,W,T + 0.5,1);
Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end
```

<div align="center"> <img src="./src/opt_pneunet.gif" width="550"> </div>


