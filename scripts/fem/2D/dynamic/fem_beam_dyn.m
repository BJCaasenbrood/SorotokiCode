%clr;
%% generate mesh from sdf
sdf = sRectangle(0,40,0,3);

msh = Mesh(sdf,'NElem',10);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/100,'TimeEnd',2,'SolverPlot',false);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);

%% select material
fem.Material = NeoHookeanMaterial(1/6,1/3); 
%% solving
fem.simulate();

%% plot energies

t = fem.Log.t;
K  = fem.Log.Kin;
U0 = fem.Log.Vg(end) + fem.Log.Psi(end);
U  = fem.Log.Psi + fem.Log.Vg - U0;
figure(103);
cla;
plot(t,K); hold on;
plot(t,U);
plot(t,U+K);
axis([0 3 0 0.02])

%% movie
close all;
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,20):numel(t)
    subplot(1,3,[1,2]);
        fem.set('Node',fem.Log.Node{ii});
        fem.show();
        axis([0 40 -40 10]);      
        
    subplot(1,3,3);
        cla;
        K = fem.Log.Kin(1:ii);
        plot(t(1:ii),K,'LineW',2); hold on;
       
        U0 = fem.Log.Vg(end) + fem.Log.Psi(end);
        U = fem.Log.Psi(1:ii) + fem.Log.Vg(1:ii) - U0;
        plot(t(1:ii),U,'LineW',2); 
        axis([0 t(end) -1e-3 1e-2]);     

        plot(t(1:ii),K+U,'LineW',2); 
    
    drawnow();
end