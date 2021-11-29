clr;
%% generate mesh
Simp  = 0.001;
GrowH = 1;
MinH  = 0.5;
MaxH  = 2;

msh = Mesh('Octarm.png','BdBox',[0,120,0,12],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.show();

%% generate fem model
fem = Fem(msh,'TimeStep',1/250,'Solver','gmres','BdBox',[0,120,-50,50]);

%% add boundary constraint
CP1 = [50,6];   % control point 1
CP2 = [120,6];  % control point 2

F1  = 0e-3;     % control force 1
F2  = 9e-3;     % control force 2

theta1 = pi/2;
theta2 = -pi/2;

fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
%fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[0,-1e-5]);
%fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Tendon',fem.FindNodes('Location',CP1),...
    F1*[cos(theta1),sin(theta1)]);

fem = fem.AddConstraint('Tendon',fem.FindNodes('Location',CP2),...
     F2*[cos(theta2),sin(theta2)]);
% 
%fem = fem.AddConstraint('Contact',@(x) SDF(x,15),[0,0]);

% assign material
fem.Material = Ecoflex0030(25);

%% solve
f = figure(101);
[~,tmpNode] = fem.solve();


%%
% t = fem.Log{1,2};
% Ux = fem.Log{2,2};
% Uy = fem.Log{3,2};
% Ve = fem.Log{12,2};
% N0 = fem.get('Node0');
% 
% figure(103); cla; subplot(1,2,1);
% for ii = 1:1:size(Ux,2)
%     Nx = N0(id,1) + Ux(:,ii);
%     Ny = N0(id,2) + Uy(:,ii)-6;
%     
%     plot(-Ny,Nx,'.','Linewidth',2,...
%         'Color',col(4,ii/(1.05*size(Ux,2))));
%     hold on;
%     Y{ii} = 1e-3*[-Ny,Nx];
% end
% 
% shp = Shapes([0,1,0,0,0,0],Y,'NModal',8);
% shp = shp.fitPOD();
% shp = shp.rebuild('NModal',1);
% 
% p = shp.string(100);
% plot(p(:,1),p(:,3)); axis equal; hold on;

%shp = shp.rebuild('NModal',1);

% 
b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'value',1,'min',1, 'max',length(tmpNode));
          
b.Callback = @(s,e) updateNode(s,e,tmpNode,[0 120 -80 50]);

function updateNode(src,~,Nd,Axis)
    class = whoClasses('Fem');
    
    low = floor(src.Value);
    upp = ceil(src.Value);

    N = lerp(Nd{low},Nd{upp},src.Value - low);
    
    for i = 1:length(class)
        class{i}.set('Node',N);
        class{i}.show('Svm');
    end
    
    axis(Axis);
end

function Dist = SDF(x,R)
Dist = dCircle(x,60,-30,R);
end
