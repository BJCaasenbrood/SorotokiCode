%% set signed distance function
sdf = @(x) Bellow(x,5,4,3,6,1.75);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,22,0,13],...
              'NElem',500,...
              'MaxIteration',50);

msh.showSDF;

%% generate mesh
msh = msh.generateMesh;
msh.show;
%% generate fem model from mesh
fem = Fem(msh);
fem = fem.set('TimeStep',1/8,...
              'ResidualNorm',1e-3);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0]);

fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,1e-3]);

fem.Material = NeoHookeanMaterial('E',5,'Nu',0.45);

%% solving
fem.solve();

function Dist = Bellow(P,x0,r1,r2,x,t)
  C1 = dCircle(P,r1+x0,0,r1);
  C2 = dCircle(P,r1+x0,0,r1-t);
  C3 = dCircle(P,r1+x0,(r1-t)*2+2*(r2),r1);
  C4 = dCircle(P,r1+x0,(r1-t)*2+2*(r2),r1-t);
  R1 = dRectangle(P,x0,x0+r1,0,r1*2+2*(r2-t));
  P1 = dIntersect(dUnion(dDiff(C3,C4),dDiff(C1,C2)),R1);
  
  R2 = dRectangle(P,r1+x0-x*0.005,r1+x0+x*1.005,r1-t,r1);
  R3 = dRectangle(P,r1+x0-x*0.005,r1+x0+x*1.005,r1+2*(r2-t),r1+2*(r2-t) + t);
  P2 = dUnion(P1,dUnion(R2,R3));
  
  C5 = dCircle(P,x0+r1+x,r1+r2-t,r2);
  C6 = dCircle(P,x0+r1+x,r1+r2-t,r2-t);
  R4 = dRectangle(P,x0+r1+x,x0+r1+x+r2,0,r1*2+2*(r2-t));
  Dist = dUnion(dIntersect(dDiff(C5,C6),R4),P2);
end
