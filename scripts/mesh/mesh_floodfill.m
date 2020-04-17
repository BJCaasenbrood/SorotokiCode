clr;
%% set signed distance function
sdf = @(x) dRectangle(x,0,1,0,1);

%% generate mesh
msh = Mesh(sdf,'BdBox',[0,1,0,1],'Quads',30);
msh = msh.generate();

%% generate density field
fem = Fem(msh);
fem.set('FilterRadius',0.01);
fem.Density = (random(0,1,fem.NElem) > 0.5)';

fem.show('E');

%% 
fem = fem.set('PressureCell',round(random(1,900,1)));
id = fem.FindElements('SDF',@(x) Subsurface(x));

msh.show('Element',id);
colormap(noir(-1))

function Dist = Subsurface(x)
Dist = dCircle(x,0.5,0.5,7.5);
end