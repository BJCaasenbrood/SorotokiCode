clr;
%% generate mesh from sdf
W = 20;
H = 1;

sdf = sRectangle(0,W,-H/2,H/2);
msh = Mesh(sdf,'NElem',100);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/150,'TimeEnd',1.5);

%% add constraint
fem = fem.addSupport(fem.FindNodes('Left'),[1,1]);
fem = fem.addGravity([0,-9.81e3]);

%% select material
fem.Material = NeoHookeanMaterial(0.1,0.3);
fem.Material.Rho = fem.Material.Rho*5;

%% solving
fem = fem.simulate();

%% extract dominated modes
NNode = 100;
Modal = [0,5,0,1,0,0];
X     = linspace(0,1,NNode)';

% shape function reconstruction
shp = Shapes(fem,Modal,'NNode',NNode,...
    'L0',W,'FilterRadius',[15,15]);

shp.g0 = SE3(eye(3),[0;0;0]);
shp = shp.reference([0,0],[W,0]);
shp = shp.reconstruct();

shp.show();
pause(2);

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;
h = [];

shp.Kp = 2;

for ii = 1:numel(t)
    
    delete(h);
    
    Nds = fem.Log.Node{ii};
    fem.set('Node',Nds);
    fem.show('Field',fem.Log.Stress{ii});
    
    Xi = shp.recoverStrain(fem,ii);
    q  = shp.estimateJointSpace(Xi,Nds);
    p  = shp.FK(q);
    
    h = fplot(p(:,[1,3]),'LineW',3,'Color',magenta);
    
    background();
    drawnow();
end
