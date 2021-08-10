clr;
%% loading .stl file
obj = Gmodel('Cylinder.stl');
obj = Blender(obj,'Center');
obj = Blender(obj,'Rotate',{'y',90});
obj = Blender(obj,'Scale',{'x',.3});

%% rendering    
obj.set('Shading','Vertex','TextureStretch',0.95);
obj.bake().render(); view(10,20);
obj.update();

axis on;
axis equal;