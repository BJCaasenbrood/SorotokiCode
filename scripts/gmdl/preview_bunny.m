clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');
obj.Texture = base;

%% rendering    
obj = obj.bake().render(); view(10,20);
obj.update();
obj.ground();
axis tight;