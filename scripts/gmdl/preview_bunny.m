clr;
%% loading .stl file
obj = Gmodel('Link3.stl');
%obj = Gmodel(@(x) dSphere(x,0,0,0,1),[-1,1,-1,1,-1,1]);
%obj.set('TextureStretch',0.625,'Shading','Face');
g = [0.9659         0    0.2588    0.0000   47.0000   11.0000  155.9064]
obj = Blender(obj,'SE3',g)
%% rendering    
obj.Texture = egg;

obj = obj.bake().render(); view(10,20);
obj.update();
axis tight;