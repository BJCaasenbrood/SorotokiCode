clr;
%% finger-actuator
T  = 1;
W1 = 5.5;
W2 = 8;
W3 = 4;
H1 = 5;
H2 = 2;
H3 = H1 - T;

R1 = sRectangle(0,W1,0,H1);
R2 = sRectangle(0,W2,-H2 - eps,0);
R3 = sRectangle(0,W3,0,H3);

S = R1 + R2 - R3;

msh = Mesh(S,'NElem',750);
msh = msh.generate();

msh.show();
%% Fem model
fem = Fem(msh,'Nonlinear',true,'TimeStep',1/35,'Linestyle','none');

%% set boundary constrains
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[0,1e-3]);

% id = fem.FindEdges('Hole',2);
% fem = fem.AddConstraint('Pressure',id,[3e-3,0]);
fem = fem.AddConstraint('Output',fem.FindNodes('Top'),[1,1]);
%% set material
fem.Material = Dragonskin10(5);

%% solving
fem.solve();

%% post-process
id = fem.FindNodes('Top');
V0 = fem.Node(id,:);

[~,I] = sort(V0(:,1));
X = V0(I,1);
Y = V0(I,2) - min(V0(:,2));

figure(103);
plot(X,Y,'-o','linewidth',2);
axis equal; axis([-6,6,0,2]);
box on; set(gca,'linewidth',1.5);
grid on;
