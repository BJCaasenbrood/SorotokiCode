clr;
%% pressure settings
P0  = 10*kpa;

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
fem = fem.set('TimeStep',1/75,'BdBox',[0,120,-80,20],'Linestyle','none',...
    'MovieAxis',[-25 120 -60 130],'Movie',0,'TimeEnd',5);

%% add boundary constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Box',[0,0,0,10]),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Pressure',fem.FindEdges('AllHole'),...
    @(t) P0*(0.4 + 0.6*sin(3*t))*min(2*t,1));

%% add output nodes
id  = fem.FindNodes('Bottom');
fem = fem.AddConstraint('Output',id,[0,0]);

%% assign material
fem.Material = Dragonskin30(25);

%% solve
out = input('* Simulate again? (y/n) ','s');
if strcmp(out,'y')
    [~,tmpNode] = fem.simulate(); 
    close all;
    save('job.mat');
else
    load('job.mat'); 
end

%% shape function reconstruction
clf;
f = figure(101);


%% post-processing
figure(105);
clf;
t   = fem.Log.t;
Nx  = fem.Log.Nx;
Ny  = fem.Log.Ny;
R   = fem.Log.Rot;
Q   = fem.Log.Str;

N      = numel(id);
ds     = 1e-3*mean(diff(Nx(1,:)));
Theta0 = [];
Gamma0 = [];

for ii = 1:numel(t)
    Rot = R{ii};
    Str = Q{ii};
    
    for jj = 1:N
        S = logm(Rot{jj});
        U = Str{jj};
        Theta0(jj,ii) = S(1,2);
        Gamma0(jj,ii) = U(1,1);
    end
    
    Theta0(:,ii) = flipud(unwrap(Theta0(:,ii)));
    Theta0(:,ii) = Theta0(:,ii) - Theta0(1,ii);
    subplot(3,1,1);
%    plot(Theta0(:,ii),...
%    'Color',col(1,min((ii +1 )/(1.00*size(Ux,2)),1))); 
%hold on;
end

%cla;
Theta = [];
Gamma = [];
for ii = 1:numel(t)    
    Gamma(:,ii) = Gamma0(:,ii);
    Theta(:,ii) = gradient(Theta0(:,ii))/ds;
%     plot(Theta(:,ii),...
%     'Color',col(1,min((ii +1 )/(1.00*size(Ux,2)),1))); 
    hold on;
end

SNAPSHOT_Q = Gamma;
SNAPSHOT_R = Theta;
[Ur,Sr,Vr]  = svd(SNAPSHOT_R);
[Uq,Sq,Vq]  = svd(SNAPSHOT_Q);

for ii = 1:2
    subplot(3,2,3);
    plot(Ur(:,ii)); hold on;
end

for ii = 1:2
    subplot(3,2,5);
    plot(Uq(:,ii)); hold on;
end

subplot(3,2,4);
plot(diag(Sr.^.5),'x-');

subplot(3,2,6);
plot(diag(Sq.^.5),'x-');

% 1;
% % figure(103); cla; 
% % subplot(1,2,1);
% % 
% % for ii = 1:1:size(Ux,2)
% % %     Nx = N0(id,1);
% % %     Ny = fem.Log.Fth(:,ii);
% %     plot(Nx,Ny,'Linewidth',3,...
% %         'Color',col(1,min((ii +1 )/(1.00*size(Ux,2)),1)));
% %     hold on;
% % end
% 
% function S = logmap(R)
% tr = trace(R);
% th = acos((tr - 1)/2);
% S = th/(2*sin(th))*(R - R.');
% end

