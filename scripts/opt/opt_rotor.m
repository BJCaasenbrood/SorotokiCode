clr;
%% generate mesh from sdf
R = 30;  % outer radius
r = 10;  % inner radius

sdf = @(x) Rotor(x,R,r);

msh = Mesh(sdf,'BdBox',[0,R,0,R],'NElem',800);
msh = msh.generate().show();

%% show generated mesh
fem = Fem(msh,'VolumeInfill',0.3,'Penal',4,'FilterRadius',R/15,...
              'Nonlinear',false,'TimeStep',1/3,...
              'OptimizationProblem','Compliant',...
              'MaxIterationMMA',50,'ChangeMax',0.15);

%% add boundary condition
id = fem.FindEdges('EdgeSelect',[0.866*R,0.5*R],20); 
fem = fem.AddConstraint('Support',id(:),[1,1]);


function D = Rotor(x,R,r)
C1 = dCircle(x,0,0,R);
C2 = dCircle(x,0,0,r);
R1 = dRectangle(x,0,R,0,R);

D = dIntersect(dDiff(C1,C2),R1);
end
