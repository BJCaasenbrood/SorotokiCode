clr;
msh = preset.mesh.diamond_bot();
msh.show('LineStyle','none');

% hold on;
% [p0, Kappa] = preset.shapes.diamond_bot_curve0();
% fplot(p0,'LineW',5,'Color',col(2));

% Kappa = Kappa(1:50);
% Kappa = [Kappa; flipud(Kappa)];

% % Kappa(isinf(Kappa)) = 0;
M = 20;
Y = pccspace(500,M);
shp = Shapes(Y,[0,M,0,0,0,0],'Length',375);
shp = shp.setBase(SE3(eye(3),[95,0,7]));
% shp.beamsolver.Xi0 = zeros(100,6);
% shp.beamsolver.Xi0(:,4) = 1;
% shp.beamsolver.Xi0(:,2) = -Kappa*0.989;

shp = shp.rebuild;
q = zeros(M,1);
q(1) = 6;
q(5) = 20;
q(6) = 5.5;
q(10) = 6;
q(11) = 6;
q(15) = 5.5;
q(16) = 20;
q(20) = 6;
G = shp.string(- q(:) * 1e-3);
p = backbone(G);
fplot(p(:,[1,3]),'LineW',3,'Color','m');