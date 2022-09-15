clr;
%% loading .stl file
obj = Gmodel('Rubberducky.stl');
obj.Texture = diffuse(0.15);

%% rendering    
obj = obj.bake().render(); view(10,20);
obj.update();
obj.ground();
axis tight;