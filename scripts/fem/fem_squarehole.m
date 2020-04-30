clr;
%% generate mesh from sdf
sdf = @(x) SquareHole(x,2,.5);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,5],'NElem',250);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'ResidualNorm',1e-4,'DisplaceNorm',1e-4,...
              'PrescribedDisplacement',true,...
              'Linestyle','none',...
              'Colormap',inferno(10));

%% add boundary condition
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,1]);

%% assign material
fem.Material = Dragonskin10A(10);

%% solving
fem.solve();

function D = SquareHole(P,H,R)
R1 = dRectangle(P,0,H/2,0,H/2);
C1 = dCircle(P,0,0,R);
D = dDiff(R1,C1);
end