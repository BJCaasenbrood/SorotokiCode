clr;
%% running fem simulation
fem = preset.fem.pneunet('contact',true);
fem.solver.TimeStep = 1/120;
fem.solver.TimeHorizon = 1;

fem = fem.simulate;

%% preset build of shapes
Y = preset.basis.chebyshev;
shp = Shapes(Y,[0,M,0,0,0,0],'Length',120);
shp = shp.setBase(SE3(eye(3),[0;0;1]));

shp = shp.reconstruct(fem);
%% plot the bending shapes
clf; I = eye(M);

for ii = 1:M
    Q = -4e-2 * I(:,ii);
    g = shp.string(Q);

    hold on; fplot(backbone(g),'LineW',3);
    view(0,0); axis equal; axis tight;
end
