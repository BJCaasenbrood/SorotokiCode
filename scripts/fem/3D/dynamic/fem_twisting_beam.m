clr;
%% generate mesh from sdf
W = 1;
H = 60;

sdf = @(x) dCube(x,-W,W,-W,W,0,H);
msh = Mesh(sdf,'BdBox',[-W,W,-W,W,0,H],'Hexahedron',[W,W,H]);

msh = msh.generate();
msh = msh.show();

%% generate fem model from mesh
fem = Fem(msh,'TimeEnd',3,'TimeStep',1/750);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1,1]);
%fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[0,0,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3,0]);

%% select material
fem.Material = Ecoflex0030(35);

%% solving
fem.simulate();

%% 
close all;
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,60):numel(t)
   % subplot(1,3,[1,2]);
        fem.set('Node',fem.Log.Node{ii});
        fem.show();
        axis([-10 10 -60 3 -30 50]);      
        view(30,10)
        box on;
    
%     subplot(1,3,3);
%         cla;
%         K  = fem.Log.Kin(1:ii);
%         U0 = fem.Log.Vg(end) + fem.Log.Psi(end);
%         U  = fem.Log.Psi(1:ii) + fem.Log.Vg(1:ii) - U0;
%         plot(t(1:ii),K,'LineW',2); hold on;
%         plot(t(1:ii),U,'LineW',2);
%         plot(t(1:ii),U+K,'LineW',2);
%         %axis([-10 10 -60 3 -30 50]);      
%         box on;    
%         
%     drawnow();

end