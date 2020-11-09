clr;
%% set signed distance function
R = 6;
sdf = @(x) dRectangle(x,-10,10,0,30);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-10,10,0,30],'Quads',25^2);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'Linestyle','none');

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Contact',@(x) SDF(x,R),[0,-0.5*R]);

%% assign material
fem.Material = Dragonskin10A;

%% solving
fem.solve();

function Dist = SDF(x,R)
Dist = dCircle(x,0,30+R,R);
end