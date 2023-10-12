clr;
%% running fem simulation
fem = preset.fem.pneunet;
fem.solver.TimeStep = 1/60;
fem.solver.TimeHorizon = 1;

fem = fem.simulate;

% preset build of shapes
M = 2;
Y = chebyspace(100,M);
shp = Shapes(Y,[0,M,0,0,0,0],'Length',120);
shp = shp.setBase(SE3(eye(3),[0;0;1]));

%% 
shp = shp.reconstruct(fem);
%% 
clf; I = eye(M);

for ii = 1:M
    Q = -4e-2 * I(:,ii);
    g = shp.string(Q);

    hold on; fplot(backbone(g),'LineW',3);
    view(0,0); axis equal; axis tight;
end
