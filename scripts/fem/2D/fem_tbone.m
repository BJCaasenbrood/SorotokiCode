clc;  clear; close all;

%% generate mesh from sdf
sdf = @(x) TensileBone(x,10,2,3,1,1);
%% 
msh = Mesh(sdf,'BdBox',[0,10,0,10],'NElem',500);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/25,'PrescribedDisplacement',true,...
         'Linestyle','none');

%% add boundary conditions
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,9]);
fem = fem.AddConstraint('Output',fem.FindNodes('Location',[1,4]),[0,1]);

%% assign material
fem.Material = Ecoflex0050();

%% solving
fem.solve();

%% plotting
% figure(101);
% subplot(2,1,2); fem.show('Svm'); view(90,90); axis tight;
% subplot(2,1,1); plot((fem.Log{2,2})/10,fem.Log{9,2},'linewidth',2,...
%     'Color',col(1)); axis tight; 
% xlabel('Uni-axial strain (-)','interpreter','latex','fontsize',12);
% ylabel('Von mises (MPa)','interpreter','latex','fontsize',12);
% grid on; set(gca,'linewidth',1);

function D = TensileBone(P,H,W,T,D,R)
dD = 0.5*(W-D);
dT = 0.5*(H-T);

R1 = dRectangle(P,0,W,0,H);
R2 = dRectangle(P,0,dD,dT,dT+T);
R3 = dRectangle(P,W-dD,W,dT,dT+T);
C1 = dCircle(P,dD-R,dT,R);
C2 = dCircle(P,dD-R,dT+T,R);
C3 = dCircle(P,W-dD+R,dT,R);
C4 = dCircle(P,W-dD+R,dT+T,R);
D0 = dDiff(dDiff(dDiff(R1,R2),C1),C2);
D = dDiff(dDiff(dDiff(D0,R3),C3),C4);
end