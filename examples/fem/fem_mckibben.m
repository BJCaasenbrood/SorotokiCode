clr;
s = sRectangle(-5,5,0,40);
msh = Mesh(s,'Quads',[10,30]);
msh = msh.generate();

E = msh.findElements('Box',[-3,3,4,36]);

fem = Fem(msh,'TimeStep',1/150);
fem = fem.addMaterial(NeoHookean);
fem = fem.addMaterial(NeoHookean(0.01,0.45));
fem = fem.setMaterial(E,2);
fem = fem.addSupport('top',[1,1]);
fem = fem.addDilation(E,0.49);

fem = fem.solve;
