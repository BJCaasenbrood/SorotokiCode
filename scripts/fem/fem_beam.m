clr;
%% generate mesh from sdf
sdf = @(x) dRectangle(x,0,5,0,1);

msh = Mesh(sdf);
msh = msh.set('BdBox',[0,5,0,1],'NElem',500);
      
msh = msh.generateMesh;

%% add boundary conditions 
fem = Fem(msh);
fem = fem.set('TimeStep',1/15);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[0,-1e-3]);
fem = fem.AddConstraint('Output',fem.FindNodes('Bottom'),[0,1e-3]);

%% select material
fem.Material = Ecoflex0030;

%% solving
fem.solve();
