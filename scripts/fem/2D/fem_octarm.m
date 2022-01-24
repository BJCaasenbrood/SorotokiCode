%clr;
%load('matlab.mat');
%% generate mesh
Simp  = 0.001;
GrowH = 1;
MinH  = 0.75;
MaxH  = 5;

msh = Mesh('Octarm.png','BdBox',[0,120,0,12],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.show();

%% generate fem model
fem = Fem(msh,'TimeStep',1/500,'BdBox',[0,120,-50,50],'TimeEnd',1.5,...
    'Linestyle','none','SolverPlot',1,'Solver','gmres');

%% add boundary constraint
CP1 = [50,6];   % control point 1
CP2 = [120,6];  % control point 2

F1  = 0e-3;     % control force 1
F2  = 5e-3;     % control force 2

theta1 = pi/2;
theta2 = -pi/2;

fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
%fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[0,-6e-4]);

%fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);

fem = fem.AddConstraint('Contact',@(x) SDF(x,10),[0,0]);

fem = fem.AddConstraint('Tendon',fem.FindNodes('Location',CP1),...
    F1*[cos(theta1),sin(theta1)]);

fem = fem.AddConstraint('Tendon',fem.FindNodes('Location',CP2),...
    F2*[cos(theta2),sin(theta2)]);

% assign material
fem.Material = NeoHookeanMaterial;

%% solve
f = figure(101);
[~,tmpNode] = fem.solve();

%%
figure(102); clf;
plot(fem.Log.t, -fem.Log.Vg,'Linewidth',4,...
              'Color',col(1)); 
hold on

H = (fem.Log.Kin);
plot(fem.Log.t,H,'Linewidth',4,...
        'Color',col(2)); 
    
plot(fem.Log.t, fem.Log.Psi,'Linewidth',4,...
              'Color',col(3)); 
hold on

H = (fem.Log.Kin) - fem.Log.Vf + fem.Log.Psi;
plot(fem.Log.t,H,'Linewidth',4,...
        'Color',col(4)); 
% %     
%hold on;
%plot(fem.Log.t,-fem.Log.Vg,'Linewidth',4,...
%           'Color',col(3)); 
%       
%plot(fem.Log.t,fem.Log.Kin,'Linewidth',4,...
%          'Color',col(4)); 
%     
xaxis('Time','-');
yaxis('Potential and Kinetic energy','J');
set(gca,'linewidth',1.5)
grid on;

%% movie
t = fem.Log.t; close all;
figure(101);


for ii = 1:fps(t,60):numel(t)
    N = tmpNode{ii};
    fem.set('Node',N);
    fem.show('Un');
    %caxis([0 1e-5]);
    axis([-60 130 -120 30]);
    drawnow();
end

function Dist = SDF(x,R)
    Dist = dCircle(x,70,-10,R);
end
