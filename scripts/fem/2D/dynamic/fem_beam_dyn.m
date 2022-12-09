%% generate mesh from sdf
sdf = sRectangle(0,40,0,3);

msh = Mesh(sdf,'Quads',[80,2]);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/250,'TimeEnd',1);
fem = fem.set('MovieAxis',[-25,40,-40,4]);

%% add constraint
fem = fem.addSupport('Left',[1,1]);
fem = fem.addGravity();
fem = fem.addContact(sCircle(5,-15,6));

%% select material
fem.Material = NeoHookeanMaterial(0.01,0.33);
fem.Material.Rr = 5;

%% solving
fig(101,[9.5,9.25]);
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
ylim([0 0.01]);

%% movie
close all;
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,300):numel(t)
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
        ylim([0 0.01]); 

        plot(t(1:ii),K+U,'LineW',2); 
    
    drawnow();
end