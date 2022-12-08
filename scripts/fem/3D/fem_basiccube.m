clr;
%% meshing
V = [-1  1  1 -1 -1  1  1 -1;
     -1 -1  1  1 -1 -1  1  1;
     -1 -1 -1 -1  1  1  1  1].';

F = [1,2,3,4,5,6,7,8];

msh = Mesh(V,F);
 
%% fem
fem = Fem(msh,'TimeStep',1/60,'Linestyle','-');
fem.Material = MooneyMaterial('C10',80,'C01',20,'K',1e7);

%%
fem = fem.addSupport(fem.FindNodes('Location',[-1,-1,-1]),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Bottom'),[0,0,1]);
fem = fem.addSupport(fem.FindNodes('Left'),[1,0,0]);
fem = fem.addDisplace(fem.FindNodes('Top'),[0,0,10]);

%%
fem = fem.solve();

%%
fig(102);
Lam = linspace(1,6,61);
Svm = horzcat(fem.Log.Stress{:});

plot(Lam,Svm(end,:))