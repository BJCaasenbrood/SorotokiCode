clr;
%% set signed distance function
sdf = @(x) dRectangle(x,-10,10,0,10);

%% generate mesh
msh = Mesh(sdf,'NElem',500,'BdBox',[-10,10,0,10]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/15,'Nonlinear',true,'FilterRadius',1e-2);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Contact',@(x) SDF(x),[0,-2]);

%% assign material
fem.Material = Dragonskin10A;
% 
%% create subsurface component
%fem.Density = fem.Density*0.25; 
% id = fem.FindElements('SDF',@(x) Subsurface(x));
% fem.Density(id) = .1;

%% solving
fem.solve();

function Dist = SDF(x)
R = 10;
Dist = dCircle(x,0,10+R,R);
end

function Dist = Subsurface(x)
%Dist = dRectangle(x,4,5.2,0,2.6);
Dist = dCircle(x,0,0,7.5);
end