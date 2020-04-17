clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,2,0,1);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,2,0,1],'Quads',[15 10]);
msh = msh.set('NElem',300,'Triangulate',false);
      
msh = msh.generate();

%% show mesh
msh.show();
set(gcf,'position',[10 10 200 200]);

%% fem
fem = Fem(msh);
fem = fem.set('TimeStep',1/10,...
              'ResidualNorm',1e-3,...
              'Nonlinear',true,...
              'FilterRadius',0.35,...
              'PrescribedDisplacement',true);
          
fem.Material = Ecoflex0030;
id = fem.FindNodes('Left');
fem = fem.AddConstraint('Support',id,[1,1]);
id = fem.FindNodes('Right');
fem = fem.AddConstraint('Support',id,[0,1]);
fem = fem.AddConstraint('Load',id,[1.0,0]);

id = fem.FindNodes('Top');
fem = fem.AddConstraint('Output',id,[0,0]);

fem.solve();
fem.show('Svm');


function d = Subsurface(x)

R = 0.75;

c1 = dCircle(x,1.75,2,R);
c2 = dCircle(x,4,2,R);
c3 = dCircle(x,6.25,2,R);
r1 = dRectangle(x,0,4,1.75,2.25);

d = dUnion(r1,dUnion(c3,dUnion(c1,c2)));

end