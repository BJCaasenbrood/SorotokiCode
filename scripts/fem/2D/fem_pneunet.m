clr;
%% pressure settings
P0  = 30*kpa;

%% generate mesh
Simp  = 0.02;
GrowH = 1;
MinH  = 1;
MaxH  = 2;

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
fem = fem.set('TimeStep',1/100,'BdBox',[0,120,-80,20],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',0.5);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
    @(t) P0*sigmoid(t));

%fem = fem.AddConstraint('Contact',@(x) SDF(x,20),[0,0]);

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

figure(103); cla; 
subplot(1,2,1);

for ii = 1:1:size(Ux,2)
    Nx = N0(id,1);
    Ny = fem.Log.Fth(:,ii);
    %Ny = Ny - Ny(end);
    %Ny = diff(Ny)*120;
    %Nx = N0(id,1) + Ux(:,ii);
    %Ny = N0(id,2) + Uy(:,ii);
    plot(Nx(1:end),Ny,'Linewidth',3,...
        'Color',col(1,min((ii +1 )/(1.00*size(Ux,2)),1)));
    hold on;
    %Y{ii} = 1e-3*[Nx,Ny];
end

% axis([-30 130 -20 130]); axis equal;
% xaxis('$x$-dimension','mm');
% yaxis('$y$-dimension','mm');
% set(gca,'linewidth',1.5)
% grid on;

% subplot(1,2,2);
% plot(t,Psi,'Linewidth',4,...
%         'Color',col(1));
% xaxis('Normalized loading','-');
% yaxis('Potential energy','J');
% set(gca,'linewidth',1.5)
% grid on;

%% movie
t = fem.Log.t; close all;
figure(105);

for ii = 1:fps(t,60):numel(t)
    N = tmpNode{ii};
    %subplot(1,2,1);
    fem.set('Node',N);
    fem.show('Svm');
    axis([-60 130 -120 30]);
    background(gitpage);
    
%     if ii == 1
%         gif('fem_pneunet_dynamic.gif','frame',gcf,'nodither');
%     else
%        gif; 
%     end
%     subplot(1,2,2);
%     plot(t(1:ii),fem.Log.Psi(1:ii),'Col',col(1),'LineW',3); hold on;
%     plot(t(1:ii),fem.Log.Kin(1:ii),'Col',col(2),'LineW',3); hold on;
%     axis([0 2 0 10]);
%     legend('Elastic','Kinetic');
%     drawnow;
    %caxis([0 1e-5]);

    drawnow();
end

%%

f = figure(102); clf;
plot(fem.Log.t, fem.Log.Psi,'Linewidth',4,...
              'Color',col(1)); 
hold on;
plot(fem.Log.t,(fem.Log.Kin),'Linewidth',4,...
        'Color',col(2)); 
    
plot(fem.Log.t, -fem.Log.Vg,'Linewidth',4,...
              'Color',col(3)); 

H = (fem.Log.Kin) - fem.Log.Vf + fem.Log.Psi;
plot(fem.Log.t,H,'Linewidth',4,...
        'Color',col(1)); 


% b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
%               'value',1,'min',1, 'max',length(tmpNode));
%           
% b.Callback = @(s,e) updateNode(s,e,tmpNode, [-60 130 -120 30]);
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
%     
%     axis(Axis);
% end
% 
function Dist = SDF(x,R)
Dist = dCircle(x,60,-R-1e-1,R);
end

