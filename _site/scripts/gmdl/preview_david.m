clr;
%% model
obj = Gmodel('David.stl');

%% set texture
obj.Texture = metal;
obj = obj.bake().render();
