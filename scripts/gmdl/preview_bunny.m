clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');
obj.Texture = softmath;
obj.TextureStretch = 0.75;
%% rendering    
obj = obj.bake().render(); view(10,20);
obj.update();
axis tight;