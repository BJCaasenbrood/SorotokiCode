clr;
%% pressure settings
P0  = 20*kpa;

%% generate mesh
msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],...
           'SimplifyTol',0.03,'Hmesh',[1,2,3],...
           'MatlabMeshType','quadratic');

msh = msh.generate();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/60,'BdBox',[-20,120,-80,20],...
              'Linestyle','none','TimeEnd',2);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,20]),[1,1]);
fem = fem.AddConstraint('Spring',fem.FindNodes('SE'),[0,5e-3]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);

fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
   P0);

%% assign material
fem.Material = Dragonskin10();
fem.Material.Zeta = 1.45;

%% solve
fem.solve(); 

%% extract dominated modes
NNode = 200;
Modal = [0,4,0,3,0,0];
X     = linspace(0,1,NNode)';

% shape function reconstruction
shp = Shapes(fem,Modal,'NNode',NNode,...
    'L0',120,'FilterRadius',[7,7]);

shp.g0 = SE3(eye(3),[0;0;0]);
shp = shp.reference([0,0],[120,0]);
shp = shp.reconstruct();

shp.show();

shppcc = Shapes(Fspace(X,NNode,8),[0,8,0,8,0,0],'L0',120,'Fem',fem);
shppcc.g0 = SE3(eye(3),[0;0;0]);
shppcc = shppcc.reference([0,0],[120,0]);
shppcc = shppcc.reconstruct();

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;
h = [];

shp.Kp = 1;
shp.Kd = 1e5;
shppcc.Kp = 1;
shppcc.Kd = 1e5;

Q = [];

for ii = numel(t)
    delete(h);
    
    Nds = fem.Log.Node{ii};
    fem.set('Node',Nds);
    fem.show('Field',fem.Log.Stress{ii});
    axis([-0 120 -55 25]);
    
    Xi = shp.recoverStrain(fem,ii);
    [q,XiE]   = shp.estimateJointSpace(Xi,Nds);
    [q0,XiE0] = shppcc.estimateJointSpace(Xi,Nds);
    
    p  = shp.FK(q);
    p0 = shppcc.FK(q0);
    
    %fplot(p0(:,[1,3]),'-','LineW',5,'Color',magenta*0.8);
    fplot(p(:,[1,3]),'LineW',5,'Color',magenta); hold on;
    
    background();
    drawnow();
end

figure(102);
clf
plot(X,Xi(:,1)*10,'LineW',1,'Color',col(1)); hold on; 
plot(X,XiE(:,1)*10,'--','LineW',2,'Color',magenta); 
plot(X,XiE0(:,1)*10,'-','LineW',2,'Color',col(3)); 
axis([0 1 -0.5 1.5])

%%
function Y = Fspace(X,N,M)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = pcc(X,ii,M); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

