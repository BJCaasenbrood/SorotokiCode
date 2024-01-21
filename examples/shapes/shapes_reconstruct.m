clr;
%% running fem simulation
fem = preset.fem.beam;
fem.solver.TimeStep = 1/120;
fem.solver.TimeHorizon = 0.5;

fem = fem.simulate;

% preset build of shapes
Y = chebyspace(100,3);
shp = Shapes(Y,[0,3,0,0,0,0],'length',100);
shp = shp.setBase(SE3(eye(3),[0;0;3]));

shp = shp.reconstruct(fem);
%% 
clf; I = eye(3);

for ii = 1:3
    Q = 2e-2 * I(:,ii);
    g = shp.string(Q);

    hold on; fplot(backbone(g),'LineW',3);
    view(0,0); axis equal; axis tight;
end
