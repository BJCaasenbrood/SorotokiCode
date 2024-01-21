clr;
%% make signed distance func.
sdf = sRectangle(10,10);

%% make mesh
msh = Mesh(sdf,'Quads',[9,9]);
msh = msh.generate();

%% show mesh;
% I = [11; 12; 13; 15; 16; 17;
%     31; 32; 33; 38; 39; 40; 41; 42; 43; 44; 
%     49; 50; 51;
%     65; 66; 67; 69; 70; 71];

% I = [14; 23;
%     31; 32; 33; 38; 39; 40; 41; 42; 43; 44; 
%     49; 50; 51;
%     59; 68];

I = [
    31; 32; 33; 40; 41; 42; 
    49; 50; 51;];

% msh = msh.remove(I);
% msh = msh.rebuild();

msh.show();
pause;
%%
fem = Fem(msh,'TimeStep',0.1);

fem = fem.addSupport('bottom',[1,1]);
id1 = fem.findNodes('location',[4.4, 10]);
id2 = fem.findNodes('location',[5.5, 10]);
fem = fem.addDisplace([id1,id2],[0 -1.0]);

fem = fem.addMaterial( NeoHookean(1,0.45) );
fem = fem.addMaterial( NeoHookean(0.001,0) );

fem = fem.setMaterial(I,2);

% solving
fem.solve('MaxIteration',5);


axis([-2 12 -2 12]);
background none;
%%
export_fig complex.png -dpng -r300