clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');
obj.Texture = grey;

%% rendering    
obj = obj.bake().render(); view(10,20);
obj.update();
axis tight;