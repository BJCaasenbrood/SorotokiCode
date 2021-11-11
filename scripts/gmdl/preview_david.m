clr;
%% model
obj = Gmodel('David.stl');

%% set texture
obj.Texture = porcelain;
obj = obj.bake().render();
