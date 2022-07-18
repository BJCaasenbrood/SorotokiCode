clr;
%% generate mesh from sdf
W = 20;
H = 1;

sdf = sRectangle(0,W,-H/2,H/2);
%msh = Mesh(sdf,'Quads',[70 2]);
msh = Mesh(sdf,'NElem',60,'Triangulate',true);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/240,'TimeEnd',0.5);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);

%% select material
fem.Material = NeoHookeanMaterial(0.01,0.3);

%% solving
fem = fem.simulate();

%% extract dominated modes
NNode = 200;
Modal = [0,3,0,1,0,0];
X     = linspace(0,1,NNode)';

% shape function reconstruction
shp = Shapes(fem,Modal,'NNode',NNode,...
    'L0',W,'FilterRadius',[35,75]);
 
shp = shp.reference([0,0],[W,0]);
shp = shp.rebuild();
shp = shp.reconstruct();

shp.show();

%%
clf;
N = 13;

Nds = fem.Log.Node{N};
Fld = fem.Log.Stress{N};
fem.set('Node',Nds);
fem.show('Field',Fld);

Xi = shp.recoverStrain(fem,N);
q  = shp.estimateJointSpace(Xi);
q(1) = q(1)*1.135;
p  = shp.FK(q);

hold on;
fplot(p(:,[1,3]),'LineW',3,'Color',magenta);

%% extract state info
% Q = fem.Log.U{1};
% 
% fem = fem.compute(Q);
% clf; 
% fem.set('Node',msh.Node);
% fem.show();
% 
% F  = fem.Log.EL.Kt*Q;
% fe = destack(fem.Log.dofa).*destack(F);
% 
% hold on;
% quiver(fem.Node(:,1),fem.Node(:,2),fe(:,1),fe(:,2));