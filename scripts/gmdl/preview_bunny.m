clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');
     
%% rendering    
obj.bake().render(); view(10,20);
obj.update();