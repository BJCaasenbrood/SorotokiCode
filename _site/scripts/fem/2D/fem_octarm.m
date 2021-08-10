clr;
%% generate mesh
Simp  = 0.001;
GrowH = 1;
MinH  = 0.15;
MaxH  = 1.5;

msh = Mesh('Octarm.png','BdBox',[0,120,0,12],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/30,'Linestyle','none','PrescribedDisplacement',true,...
    'Movie',true,'MovieAxis',[0 130 -50 20]);

%% add boundary constraint
CP1 = [60,10.5];  % control point 1
CP2 = [60,1.5];   % control point 1
CP3 = [120,8.7];  % control point 1
CP4 = [120,3.3];  % control point 1
F1  = 30e-3;      % control force 1
F2  = 0;          % control force 1
F4  = -10e-3;      % control force 1
F3  = 0;          % control force 1

fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',[120,8.9]),[-62,-18]);
%fem = fem.AddConstraint('Load',fem.FindNodes('Location',[120,3.3]),[-60,-15]);
% fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP1),[0,F1]);
% fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP2),[0,F2]);
% fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP3),[0,F3]);
% fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP4),[F4,0]);

%% output nodes
% id = fem.FindNodes('Box',[0 120 5.5 6.5]);
% fem = fem.AddConstraint('Output',id,[0,0]);
% 
% %% assign material
fem.Material = Dragonskin10;
% 
% %% solve
f = figure(101);
[~,tmpNode] = fem.solve();
% 
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
% b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
%               'value',1,'min',1, 'max',length(tmpNode));
%           
% b.Callback = @(s,e) updateNode(s,e,tmpNode,[0 120 -50 50]);
% 
% function updateNode(src,~,Nd,Axis)
%     class = whoClasses('Fem');
%     
%     low = floor(src.Value);
%     upp = ceil(src.Value);
% 
%     N = lerp(Nd{low},Nd{upp},src.Value - low);
%     
%     for i = 1:length(class)
%         class{i}.set('Node',N);
%         class{i}.show('Un');
%     end
%     axis(Axis);
% end

