clr;
%% loading .stl file
obj = Gmodel('Fandisk.stl');
%obj = Gmodel(@(x) dSphere(x,0,0,0,1),[-1,1,-1,1,-1,1]);
obj.set('TextureStretch',0.625,'Shading','Face');

%% rendering    
obj.Texture = bump;

obj = obj.bake().render(); view(10,20);
obj.update();