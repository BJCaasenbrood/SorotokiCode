clr;
%% set signed distance function
R = 4;
sdf = sRectangle(0,10,0,10);

%% generate mesh
msh = Mesh(sdf,'NElem',250);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/500,'Linestyle','none',...
    'PrescribedDisplacement',1);

%% add constraint
fem = fem.addSupport(fem.FindNodes('Left'),[1,0]);
fem = fem.addSupport(fem.FindNodes('Bottom'),[0,1]);
fem = fem.addContact(SDF(R),[0,-1.5*R]);

%% assign material
fem.Material = Dragonskin10(40);

%% solving
fem.solve();

function sdf = SDF(R)
sdf = sCircle(0,10+R,R);
end