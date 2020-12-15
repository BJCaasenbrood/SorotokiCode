---
layout: default
title: Topology optimization
parent: Documentation
nav_order: 3
---

# Topology optimization


# Numerical examples
### Example: Pneu-net soft robot

<div align="center"> <img src="./img/opt_pneunet.gif" width="550"> </div>

```matlab
%% generate mesh from sdf
sdf = @(x) PneuNet(x,20,40,1,20);

msh = Mesh(sdf,'BdBox',[0,20,0,40],'Quads',[25 50]);
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',4,...
              'Nonlinear',false,'TimeStep',1/3,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',70);

%% set spatial settings
fem = fem.set('Periodic',[0.5, 0],'Repeat',ones(8,1));

%% add boundary condition
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = fem.FindNodes('Right'); 
fem = fem.AddConstraint('Spring',id,[0,1]);
fem = fem.AddConstraint('Output',id,[0,-1]);

id = fem.FindElements('Location',[10,25],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[10,25],0.5);

%% material
fem.Material = Dragonskin10;

%% solving
fem.optimize();
fem.show('ISO');

function Dist = PneuNet(P,W,H,E,T)
R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,-W/2,E,T,H+H/2);
R3 = dRectangle(P,W-E,W+W/2,T,H+H/2);
C1 = dCircle(P,0,T + 0.5,1);
C2 = dCircle(P,W,T + 0.5,1);
Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end
```

### Example: Pneumatic bellow

<div align="center"> <img src="./img/opt_bellow.gif" width="350"> </div>


```matlab
%% set signed distance function
W = 8;
H = 4;
sdf = @(x) Bellow(x,W,H);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,W,0,H],'NElem',750);
msh = msh.generate();

%% generate fem from mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',0.75,...
              'Nonlinear',false,'TimeStep',1/3,'ReflectionPlane',[1,1],...
              'OptimizationProblem','Compliant','Repeat',[1 2],...
              'MaxIterationMMA',65);

%% add boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);

id = fem.FindNodes('Location',[0.01*W,H]);
fem = fem.AddConstraint('Output',id,[0,1]);
fem = fem.AddConstraint('Spring',id,[0,.1]);

id = fem.FindNodes('Line',[0.02*W,W,H,H]);
fem = fem.AddConstraint('Spring',id,[0,.1]*1e-1);

id = fem.FindElements('Location',[0,0],1);
fem = fem.AddConstraint('PressureCell',id,[1e-3,0]);

%% set density
fem = fem.initialTopology('Hole',[0,0],1.0);

%% material
fem.Material = Ecoflex0030;

%% solving
fem.optimize();

function D = Bellow(x,W,H)
R1 = dRectangle(x,0,W,0,H);
C2 = dCircle(x,W,H,1.);

D = dDiff(R1,C2);
end
```



