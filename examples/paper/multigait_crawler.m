clr;

% load preset fem
fem = preset.fem.multigait_crawler('n',1250);

fem = fem.addPressure(fem.findEdges('BoxHole',[0 50 0 15]),...
   @(t) pset(t,1));
fem = fem.addPressure(fem.findEdges('BoxHole',[50 100 0 15]),...
   @(t) 0.75*pset(t,2));
fem = fem.addPressure(fem.findEdges('BoxHole',[100 150 0 15]),...
   @(t) pset(t,3));

% solve dynamic simulation
fem = fem.simulate();

function y = pset(t,id)
    w  = 2*pi;
    P0 = 15 * 1e-3;
    t  = mod(t,.75);
    y  = P0*tsin(w*t - (id-1)*pi/5).*sigmoid(3*w*t);
end
