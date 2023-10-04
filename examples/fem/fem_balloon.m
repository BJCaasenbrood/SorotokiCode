clr;
T = 1;
R = 2;
L = 20;
sdf = sRectangle(L,2 * R) + sCircle(R,[L,R]) ... 
    - sRectangle([T,T],[L,T+R]) - sCircle(R-T,[L,R]);

msh = Mesh(sdf,'NElem',175);
msh = msh.generate();

fem = Fem(msh,'LineStyle','none');
fem = fem.addMaterial( NeoHookean(0.01,0.33) );

fem = fem.addSupport('left',[1,1]);
fem = fem.addPressure('allhole',1e-3);

fem.solve('TimeStep',1/3333);