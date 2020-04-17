clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,8,0,4);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,8,0,4]);
msh = msh.set('NElem',2000,'Triangulate',false);
      
msh = msh.generate();

%% show mesh
msh.show();
set(gcf,'position',[10 10 500 500]);
colormap(bluesea);
caxis([-4 1]);

%% fem
fem = Fem(msh);
fem = fem.set('TimeStep',1/10,...
              'ResidualNorm',1e-3,...
              'Nonlinear',false,...
              'FilterRadius',0.35,...
              'PrescribedDisplacement',true);
          
% fem.Material = Ecoflex0030;
% id = fem.FindNodes('SW');
% fem = fem.AddConstraint('Support',id,[1,1]);
% id = fem.FindNodes('SE');
% fem = fem.AddConstraint('Support',id,[1,1]);
% 
% id = fem.FindNodes('Location',[1,1]);
% fem = fem.AddConstraint('Load',id,[0,-0.25]);
% 
% id = fem.FindNodes('Top');
% fem = fem.AddConstraint('Output',id,[0,0]);
% 
% fem.optimize();
% fem.show('E')

id = fem.FindElements('SDF',@(x) Subsurface(x));
fem.Density(id) = .01;
caxis([0 1])

fem.show('E');

id = fem.FindElements('Location',[0,2]);
fem = fem.AddConstraint('PressureCell',id,[0,0]);

id = fem.FindElements('FloodFill',fem,fem.Density);

hold on
msh.show('Element',id);
colormap(metro(-1))


function d = Subsurface(x)

R = 0.75;

c1 = dCircle(x,1.75,2,R);
c2 = dCircle(x,4,2,R);
c3 = dCircle(x,6.25,2,R);
r1 = dRectangle(x,0,4,1.75,2.25);

d = dUnion(r1,dUnion(c3,dUnion(c1,c2)));

end