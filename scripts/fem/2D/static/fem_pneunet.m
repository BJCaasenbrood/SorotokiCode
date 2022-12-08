clr;
%% simulation settings
L = 120;  % Lenght
H = 20;   % Height
T = 15;   % Thickness
P = 5 * kpa;

%% finite element settings
Simp  = 0.02;
GrowH = 1;
MinH  = 2;
MaxH  = 3;

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,L,0,H],...
    'SimplifyTol',Simp,'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.generate();

%% generate fem model
fem = Fem(msh,'TimeStep',1/1000);

%% add boundary constraint
fem = fem.addSupport('Left',[1,1]);
fem = fem.addGravity();
fem = fem.addPressure([],P);

%% assign material
fem.Material = Ecoflex0030;%Dragonskin10();

%% solve
fem.solve();

%% movie
figure(105);
t = fem.Log.t; 

Pos = [];

for ii = 1:numel(t)
    Pos(:,:,ii) = fem.Log.Node{ii}(fem.FindNodes('Bottom'),:);
    
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    %axis([-30 130 -100 30]);   
    background();
    drawnow();
end

%%
figure(106); clf;
h = msh.show();
%set(h,'Linestyle','none');

C = cmap_orange(numel(t));
for ii = 1:numel(t)
    plot(Pos(:,1,ii),Pos(:,2,ii),...
        'LineW',2,'Color',C(ii,:)); hold on;
end

plot(reshape(Pos(1,1,:),[217 1]),reshape(Pos(1,2,:),[217 1]),...
        'LineW',2,'Color',C(end,:)); hold on;
    
axis equal; axis([-30 130 -100 30]);   
axis on;

