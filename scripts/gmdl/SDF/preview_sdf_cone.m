clr;
%% model
Dist = @(X) dCone(X,0,0,0,1,2);
obj = Gmodel(Dist,domain(-2,2,3),'Quality',120,'Shading','Face');
%% set texture
obj.Texture = base;
obj.bake();

%% show
obj.render(); 