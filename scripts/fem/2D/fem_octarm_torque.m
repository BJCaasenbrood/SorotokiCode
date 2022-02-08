clr;
%% generate mesh
Simp  = 0.001;
GrowH = 1;
MinH  = 1;
MaxH  = 2;

msh = Mesh('Octarm.png','BdBox',[0,120,0,12],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/10);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Displace',fem.FindNodes('Right'),so2(0.5*pi));

% assign material
fem.Material = NeoHookeanMaterial;
%Dragonskin30(25);

%% solve
f = figure(101);
[~,tmpNode] = fem.solve();

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
          
b.Callback = @(s,e) updateNode(s,e,tmpNode,[0 120 -50 50]);

function updateNode(src,~,Nd,Axis)
    class = whoClasses('Fem');
    
    low = floor(src.Value);
    upp = ceil(src.Value);

    N = lerp(Nd{low},Nd{upp},src.Value - low);
    
    for i = 1:length(class)
        class{i}.set('Node',N);
        class{i}.show('Un');
    end
    
    axis(Axis);
end

