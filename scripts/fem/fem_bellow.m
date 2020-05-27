clr;
%% set signed distance function
sdf = @(x) Bellow(x,5,4,6,5,7,5,2);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,25,0,25],'NElem',500);
msh = msh.generate();

%% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/50,'PrescribedDisplacement',true,...
    'SigmoidFactor',0.5);

%% add constraint
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Top'),[1,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,5]);

%% add logger nodes
fem = fem.AddConstraint('Output',fem.FindNodes('Location',[6,22]),[0,0]);

%% assign material
fem.Material = Dragonskin10A;

%% solving
fem.solve();

%% plot nonlinear stiffness
figure(101);
uy = mean(fem.Log{3,2},2);
fy = mean(fem.Log{6,2},2)*1e3;

subplot(1,2,1); fem.show();
subplot(1,2,2); plot(uy,diff([0;fy])./diff([0;uy]),'linewidth',2,'Color',col(2));
xlabel('Displacement (mm)','interpreter','latex','fontsize',12);
ylabel('Reaction force (mN)','interpreter','latex','fontsize',12);
grid on; set(gca,'linewidth',1); axis square;

function Dist = Bellow(P,r0,r1,r2,r3,r4,x,t)
  C1 = dCircle(P,r2+r0,0,r1);
  C2 = dCircle(P,r2+r0,0,r2);
  R1 = dRectangle(P,r0,r2+r0,0,r2);
  R2 = dRectangle(P,r2+r0,r2+x+r0,r1,r2);
  C3 = dCircle(P,r2+x+r0,r2+r3,r3);
  C4 = dCircle(P,r2+x+r0,r2+r3,r4);
  R3 = dRectangle(P,r2+x+r0,r2+x+r4+r0,r1,r1+2*r4);
  R4 = dRectangle(P,r2+r0,r2+x+r0,r2+2*r3,r2+2*r3+t);
  R6 = dRectangle(P,r0,r2+r0,r2+2*r3,r2+r2+2*r3);
  C5 = dCircle(P,r2+r0,2*r2+2*r3,r2);
  C6 = dCircle(P,r2+r0,2*r2+2*r3,r1);
  
  Dist0 = dDiff(dIntersect(R6,C5),C6);
  Dist1 = dUnion(dDiff(dIntersect(R3,C4),C3),R4);
  Dist  = dUnion(dUnion(dDiff(dUnion(R2,dIntersect(R1,...
      C2)),C1),Dist1),Dist0);
end
