clr;
%% model
obj = Gmodel('SoftRoboticArm.stl');

%% set texture
obj.Texture = mateplastic;
obj = obj.bake();

%% show
obj = obj.render();
view(130,21);
