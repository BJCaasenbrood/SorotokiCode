clr;
%% 
L = 250;    % length of robot
M = 8;      % number of modes
N = 75;     % number of discrete points on curve
H = 1/50;   % timesteps
FPS = 50;   % animation speed

%% 
Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,0,0,0,0],'L0',L);
shp.g0 = SE3(roty(pi/2),zeros(3,1));

%%
mdl = Model(shp,'TimeEnd',30);
mdl.X0(1) = 0.05;
mdl = mdl.simulate(); 

%% 
figure(100);
plot(mdl.Log.t,mdl.Log.x(:,1:M),'LineW',2);
colororder(col);

%% animation
for ii = 1:fps(mdl.Log.t,FPS):length(mdl.Log.t)
    shp = shp.render(mdl.Log.x(ii,1:M));
    
    axis([-.5*L .5*L -.5*L .5*L -L 0.1*L]);
    view(30,30);
    drawnow();
end
