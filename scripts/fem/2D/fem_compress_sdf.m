clr;
%% set signed distance function
R = 6;
sdf = @(x) dRectangle(x,-10,10,0,30);

%% generate mesh
msh = Mesh(sdf,'BdBox',[-10,10,0,30],'Quads',[10,50]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/55,'Linestyle','none','Solver','lu');

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Contact',@(x) SDF(x,R),[0,-0.5*R]);

%% assign material
fem.Material = Dragonskin10;

%% solving
fem.solve();

function Dist = SDF(x,R)
Dist = dCircle(x,0,30+R,R);
end