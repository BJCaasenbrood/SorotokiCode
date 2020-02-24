clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,10,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,10,0,10],'Quads',10);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/150,'Nonlinear',true);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Contact',@(x) SDF(x),[0,-5]);

%% assign material
fem.Material = Ecoflex0030;

%% solving
fem.solve();

function Dist = SDF(x)
Dist = dCircle(x,0,15,5);
end