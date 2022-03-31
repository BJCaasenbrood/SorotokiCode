clr;
%% generate mesh
Simp  = 0.001;
GrowH = 1;
MinH  = 0.75;
MaxH  = 5;

msh = Mesh('Octarm.png','BdBox',[0,120,0,12],...
    'SimplifyTol',Simp,'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.show();

%% find tendon attachments
CP1 = [50,6];  % control point 1
CP2 = [100,6];  % control point 2

id1 = FindNode(msh.Node,'Location',CP1);
id2 = FindNode(msh.Node,'Location',CP2);

msh.show('Node',[id1,id2]); pause(2);
%% generate fem model
fem = Fem(msh,'TimeStep',1/250,'TimeEnd',5,'Linestyle','none');

%% add boundary constraint
F1  = 1e-2;     % control force 1
F2  = 1e-3;      % control force 2

theta1 = pi/2;
theta2 = -pi/2;

fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);

fem = fem.AddConstraint('Tendon',fem.FindNodes('Location',CP1),...
    F1*[cos(theta1),sin(theta1)]);

fem = fem.AddConstraint('Tendon',fem.FindNodes('Location',CP2),...
    F1*[cos(theta2),sin(theta2)]);

%% assign material
fem.Material = Dragonskin10(30);
fem.Material.Zeta = 0.75; 

%% solve
fem = fem.simulate();

%% movie
t = fem.Log.t; 
close all;
figure(101);

for ii = 1:fps(t,120):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    axis([-60 130 -120 30]);
    drawnow();
end