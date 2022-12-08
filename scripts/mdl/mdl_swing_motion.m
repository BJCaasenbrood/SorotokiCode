clr;
%%
Omega = 0.5*pi; % rad/s
Ampli = 235;    % Nmm (N-millimeter)

%% 
L   = 100;     % length of robot
M   = 6;       % number of modes
N   = 30;      % number of discrete points on curve
FPS = 150;     % animation speed

%% 
Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,0,0,0,0],'Length',L);
shp = shp.addGravity([0;0;9810]);
shp = shp.setRadius(8);
shp = shp.setRamp(.7);
shp.Material = NeoHookeanMaterial(.25,0.33);

shp = shp.rebuild();

%% simulate model
mdl = Model(shp,'TimeEnd',15,'TimeStep',1/60);

mdl.Controller = @(M) Controller(M,Omega,Ampli);
mdl = mdl.simulate(); 

%% animation
%[rig] = setupRig(shp);
t = mdl.Log.t;
for ii = 1:fps(t,FPS):length(t)

    shp = shp.render(mdl.Log.x(ii,1:M));
    
    axis([-.5*L .5*L -.5*L .5*L -0.1*L L]);
    view(30,30);
    drawnow();
end

%% setup controller
function tau = Controller(mdl,w,A)
    tau = zeros(mdl.NDim/2,1);
    tau(1) = A*sin(w*mdl.t);
    tau(2) = 0.1*A*sin(w*mdl.t);
end
