clr;
%%
obj = Gmodel('gripperbase.stl');
obj.Texture = grey*1.5;
obj.bake();
obj.render();

obj2 = Gmodel('gripperholder.stl','Shading','Face');
obj2.Texture = whitebase*0.65;
obj2.render();

sdf = sSphere(0,0,-60,45); hold on;
sdf.show('Quality',25); 
%%
%
M = 5;
W = 5;
T = 0.85;
L = 130;
[Y,X] = chebyspace(100,M);
Y = Y.*smoothstep(X) + 1e-6;

shp = Shapes(Y,[M,M,0,0,0,0],'L0',L,...
    'xia0',[0,0,0,0,0,1],'Texture',diffuse(.45)*1.15);
shp = shp.setRadius([4,1,1]);
shp.Kp = 30;
shp.Kd = 1e8;
%%

f = fig(101,[9.5,10]);
G{1} = SE3(rotx(pi + 0.175),[18,20,20]);
G{2} = SE3(rotx(pi + 0.175),[-18,20,20]);
G{3} = SE3(rotz(pi/2)*rotx(pi),[-39,0,20]);
G{4} = SE3(rotz(-pi/2)*rotx(pi),[39,0,20]);
G{5} = SE3(rotx(-pi - 0.175),[18,-20,20]);
G{6} = SE3(rotx(-pi - 0.175),[-18,-20,20]);

%Q = zeros(6,shp.NJoint);
Q      = zeros(6,shp.NJoint) ;
% Q(:,1) = -1;
% Q(:,2) = -10;
% Q(1,1) = 1;
% Q(1,2) = -15;
% Q(2,1) = 14;
% Q(2,2) = 0;
% Q(3,1) = 2;
% Q(3,3) = -8;
% Q(4,1) = 2;
% Q(4,3) = -8;
% Q(5,1) = 8;
% Q(5,2) = 8;
% Q(6,1) = 20;
% Q(6,2) = -1;
% Q(:,2) = -10*kpa;
% Q(:,3) = -10*kpa;

% Q = Q*kpa;
% 
% for ii = 1:6
%     
%     shp.g0 = G{ii};
%     view(190,-10);
%     h{ii} = shp.render(Q(ii,:)*kpa);
% end
% % % 
% axis tight;
% drawnow;
% background();
% zoom(1.35);
% for ii = 1:6
% delete(h{i});
% end

for ii = 1:6
    
    shp.g0 = G{ii};
    %shp = shp.rebuild();
       
    view(190,-10);
    V  = shp.FK(Q(ii,:));
    P  = sdf.surfaceproject(V);
    qd = shp.IK(P,Q(ii,:));

    shp.render(qd + 1e-6);
end

obj.update();
obj2.update();
axis tight;

drawnow;
background();
zoom(1.25);

% %%
% I = frame2im(getframe(gcf));
% I = imrotatepad(I,-10);
% % clf;
% % imshow(I)
% %%
% close all;
% f = fig(101,[10,9.5]); delete(gca);
% ha = tight_subplot(1,2,[.01 .06],[.01 .01],[.01 .01]);
% axes(ha(1));
% 
% im = imread('foo2.png');
% imshow(im);
% 
% axes(ha(2));
% imshow(I);
% 
% %%
% figpath = @(x) ['~/Documents/phd/papers/2023/Sorotoki/IEEEAcces/fig/',x];
% % 
% W0 = 0.400;
% X0 = num2str(W0,4);
% Y0 = num2str((f.InnerPosition(4)/f.InnerPosition(3))*W0,4);
% % 
% cleanfigure('targetResolution',100);
% matlab2tikz(figpath('fig_mdl_example2.tex'),...
%       'width',[X0,'\textwidth'],'height',[Y0,'\textwidth']);
