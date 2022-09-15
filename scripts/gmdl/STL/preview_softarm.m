clr;
%% loading .stl file
obj = Gmodel('SoftRoboticArm_.stl','Shading','Face');
obj.Texture = diffuse(0.5);

%% rendering    
obj = obj.bake().render(); view(10,20);
obj.update();
axis tight;