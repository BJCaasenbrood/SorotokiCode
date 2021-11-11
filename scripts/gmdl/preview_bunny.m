clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');

%% rendering    
obj.Texture = bump;
obj.set('TextureStretch',0.99);

obj = obj.bake().render(); view(10,20);
obj.update();
axis tight;