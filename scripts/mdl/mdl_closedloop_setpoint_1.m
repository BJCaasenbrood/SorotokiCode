clr; 
%% 
L = 120;   % length of robot
M = 10;     % number of modes
N = 120;    % number of discrete points on curve
H = 1/60;  % timesteps
FPS = 30;  % animation speed

Modes = [0,M,M,0,0,0];  % pure-XY curvature

%% shapes
% generate nodal space
Y = chebyspace(N,M);
shp = Shapes(Y,Modes,'Length',L,'RenderMuscles',0,...
    'xia0',[0,0,0,1,0,0],...
    'Texture',softmath);

%% set material properties
shp = shp.setRadius(10);
shp = shp.setRamp(0.75);
shp.Material = NeoHookeanMaterial(2,0.49);

shp.Gravity = [0;0;-9810];

Nc = 0.5;
Muscle = @(s,N) [0*s;cos(N*pi/2 + 2*pi*Nc*s); sin(N*pi/2+ 2*pi*Nc*s)];
shp = shp.addMuscle(@(s) Muscle(s,1e-3));

shp = shp.rebuild(); 

%% build model class
mdl = Model(shp,'TimeEnd',2.5,'MaxIteration',10);

%% controller
mdl.Controller = @(M) Controller(M);

%% simulate system
mdl = mdl.simulate(); 

% %%
% pause(0.1);
% thesispath = @(x) ['/home/brandon/Documents/phd/thesis/3_chapters/',x];
% 
% W0 = 0.875;
% X0 = num2str(W0,4);
% Y0 = num2str((f.InnerPosition(4)/f.InnerPosition(3))*W0,4);
% 
% cleanfigure('targetResolution', 100);
% matlab2tikz(thesispath('3_chapter/img/fig_C3_modes_converge.tex'),...
%      'width',[X0,'\textwidth'],'height',[Y0,'\textwidth']);

%% animation
sph = Gmodel(sSphere(4));
sph.Texture = diffuse(0.925);
%sph.bake.render();
fig(101,[9.5,9]);
set(gcf,'GraphicsSmoothing','off')
jj = 1;
for ii = [1,38,55,90,150];

    view(30,30);
    shp = shp.render(mdl.Log.x(ii,1:2*M));
    
    axis([-.1*L L -.1*L .5*L -0.1*L 0.1*L]);
    zoom(1.5)
    drawnow();
    
    export_fig(gca,['img',num2str(jj),'.png'],'-nocrop','-q1200',...
        '-painters','-a4')
    jj = jj + 1;
end


function gd = gref(t)
gd = SE3(eye(3),[50 + 30*sin(t),10*cos(t),10]);
end

%% setup controller
function tau = Controller(mdl)
t = mdl.t;
% q  = mdl.Systems{1}.Log.q;
% gg = mdl.Systems{1}.Log.FK.g;
% JJ = mdl.Systems{1}.Log.FK.J;
% 
% ge = gg(:,:,end);
% J  = admapinv(ge)*JJ(:,:,end);
% gd = gref(0);
% 
% lam1 = 5e4;
% lam2 = 1e4;
% Kp = diag([1e-5,1e-5,1e-5,1,1,1]);
% 
% Xi = smoothstep(4*t)*logmapSE3(ge\gd);
% Fu = Kp*tmapSE3(Xi)*wedge(Xi);
% 
% tau = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
% tau = tau + mdl.Systems{1}.Log.EL.fg + mdl.Systems{1}.Log.EL.K*q;
tau = mdl.Systems{1}.Log.EL.fg*0;
tau(1) = -sigmoid(t)*75;
end


