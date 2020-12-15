clr;
%% model
Dist = @(X) dCube(X,-1,1,-1,1,-1,1);
obj = Gmodel(Dist,domain(-2,2,3),'Quality',80,'Shading','Face');
%% set texture
obj.Texture = base;
obj.bake();

%% show
obj.render(); 