clr;
%% settings
P0  = 30*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;
MinH  = 2;
MaxH  = 3;

msh = Mesh('Pneunet.png','BdBox',[0,120,0,20],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

figure(101);
subplot(2,1,1); imshow('Pneunet.png');
subplot(2,1,2); msh.show();

%% re-orient the mesh
%msh = Blender(msh,'Rotate',{'x',90}); msh = msh.show(); pause(2);

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/200,'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',.2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),[P0,0]);

%% add output nodes
id  = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = Dragonskin30(25);

%% solve
clf;
f = figure(101);
[~,tmpNode] = fem.simulate();

%% post-processing
t   = fem.Log.t;
Ux  = fem.Log.Ux;
Uy  = fem.Log.Uy;
Psi = fem.Log.Psi;
N0  = fem.get('Node0');

figure(103); cla; subplot(1,2,1);
for ii = 1:1:size(Ux,2)
    Nx = N0(id,1) + Ux(:,ii);
    Ny = N0(id,2) + Uy(:,ii);
    plot(Nx,Ny,'Linewidth',3,...
        'Color',col(4,(ii)/(1.00*size(Ux,2))));
    hold on;
    Y{ii} = 1e-3*[Nx,Ny];
end

axis([-30 130 -20 130]); axis equal;
xaxis('$x$-dimension','mm');
yaxis('$y$-dimension','mm');
set(gca,'linewidth',1.5)
grid on;

subplot(1,2,2);
plot(t,Psi,'Linewidth',4,...
        'Color',col(1));
xaxis('Normalized loading','-');
yaxis('Potential energy','J');
set(gca,'linewidth',1.5)
grid on;

%% movie
t = fem.Log.t; close all;
figure(101);

for ii = 1:fps(t,200):numel(t)
    N = tmpNode{ii};
    fem.set('Node',N);
    fem.show('Un');
    %caxis([0 1e-5]);
    axis([-60 130 -120 30]);
    drawnow();
end



% b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
%               'value',1,'min',1, 'max',length(tmpNode));
%           
% b.Callback = @(s,e) updateNode(s,e,tmpNode,[0 140 -120 120]);
% 
% function updateNode(src,~,N
% 
% d,Axis)
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
%     
%     axis(Axis);
% end
% 
