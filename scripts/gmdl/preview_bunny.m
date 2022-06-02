clr;
%% loading .stl file
obj = Gmodel('Bunny.stl','Alpha',0.2);
obj.Texture = bubble;

%% rendering    
obj = obj.bake().render(); view(10,20);
obj.update();
axis tight;