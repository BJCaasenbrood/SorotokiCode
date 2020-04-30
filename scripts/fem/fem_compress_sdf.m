clr;
%% set signed distance function
sdf = @(x) dRectangle(x,-10,10,0,10);

%% generate mesh
msh = Mesh(sdf,'NElem',250,'BdBox',[-10,10,0,10]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'Nonlinear',true,'FilterRadius',1e-2);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Contact',@(x) SDF(x),[0,-2]);

%% assign material
fem.Material = Dragonskin10A;

%% solving
fem.solve();

function Dist = SDF(x)
R = 3;
Dist = dCircle(x,0,10+R,R);
end