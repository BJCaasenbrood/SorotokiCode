clr;
%% model
obj = Gmodel('Fandisk.stl');
obj.set('Shading','Face');
%% set texture
obj.Texture = diffuse(0.95);
obj.Emission = [0.15 0.75 0.15];

%% set view
obj.bake().render().update();
obj.ground();
axis tight;
view(30,30);