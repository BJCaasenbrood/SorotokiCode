<div align="center"> <img src="./src/fem.png" width="650"> </div>

# Nonlinear Topology Optimization

[**Homepage**](https://bjcaasenbrood.github.io/SorotokiCode/)

### Example: Pneu-net soft robot
```matlab
%% set signed distance function
sdf = @(x) PneuNet(x);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,20,0,20],...
              'NElem',500);
      
msh = msh.generateMesh;

%% show generated mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/5,...
              'VolumeInfill',0.3,...
              'VolumetricPressure',true,...
              'FilterRadius',1.5,...
              'Periodic',[20, 0],...
              'Repeat',[1, 1, 1, 1],...
              'OptimizationProblem','Compliant');

%% add constraint
id = fem.FindNodes('Left'); 
fem = fem.AddConstraint('Support',id,[1,1]);

id = [fem.FindNodes('Right'); fem.FindNodes('Line',[15,20,0,0])]; 
fem = fem.AddConstraint('Output',id,[0,-1]);
fem = fem.AddConstraint('Spring',id,[0,2]);

id = fem.FindElements('Location',[10,10],1);
fem = fem.AddConstraint('PressureCell',id,[1e-2,0]);

%% set density
fem = fem.initialTopology('Hole',[10,10],2);

%% material
fem.Material = Dragonskin10A;

%% solving
fem.optimize();

%% signed distance function for Pneu-net
function Dist = PneuNet(P)
  R1 = dRectangle(P,0,20,0,20);
  R2 = dRectangle(P,-4,1,4,20);
  R3 = dRectangle(P,19,24,4,20);
  C1 = dCircle(P,0,4.5,1);
  C2 = dCircle(P,20,4.5,1);
  Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end


```

