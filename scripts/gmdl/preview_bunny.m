clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');

%% rendering    
fig(101,[9,9]);
obj = obj.bake().render(); 
axis tight;

%% material
obj.Texture = bump;
obj.update();