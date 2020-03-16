clr;
%% set signed distance function
sdf = @(x) dRectangle(x,-10,10,0,10);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-10,10,0,10],'Quads',50);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'Nonlinear',true,'FilterRadius',1e-2);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Contact',@(x) SDF(x),[0,-1]);

%% assign material
fem.Material = Dragonskin10A;

%% create subsurface component
fem.Density = fem.Density*0.25; 
id = fem.FindElements('SDF',@(x) Subsurface(x));
fem.Density(id) = 1;

%% solving
fem.solve();

function Dist = SDF(x)
Dist = dCircle(x,0,15,5);
end

function Dist = Subsurface(x)
Dist = dRectangle(x,4,5.2,0,2.6);
end