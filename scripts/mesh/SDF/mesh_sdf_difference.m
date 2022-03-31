clr;
%% generate SDF class
sdf1 = sCircle(0, 0, 0.5);
sdf2 = sCircle(0.5, 0, 0.5);

%% math operations
f = sdf1-sdf2;
f.cmap = turbo;
f.show(); view(0,90);
pause(1);

%% meshing
msh = Mesh(f);
msh = msh.generate();
msh.show();

figure(101);
subplot(1,2,1); f.show(); 
subplot(1,2,2); msh.show();