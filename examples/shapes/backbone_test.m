clr;
% set parameters
[N,M] = deal(54,2);

% build polyspace
Y = chebyspace(N,M);

% construct Shapes class
shp = Shapes(Y,[0,M,M,0,0,0],'Quality',175);
shp = shp.setRamp(0.8);  % set ramp

Q = random(-1,1,M*2) * 6e-2;

g = shp.string(Q);
SE3plot(g);

view(30,30);
axis equal; axis tight;
box on; 

function SE3plot(g)
    [p,ux,uy,uz] = backbone(g);
    fplot(p,'-','LineW',1); hold on;
    fquiver(p,ux,0.1,'Color',col(1),'LineW',2);
    fquiver(p,uy,0.1,'Color',col(2),'LineW',2);
    fquiver(p,uz,0.1,'Color',col(3),'LineW',2);
end