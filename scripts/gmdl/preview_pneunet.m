clr;
%% model
obj = Gmodel('Pneunet.stl');
obj.set('TextureStretch',0.95);

%% deform
Blender(obj,'Translate',{'y',-0.4});
Blender(obj,'Curve',{'PCC+',120,0,0.8});
Blender(obj,'Rotate',{'y',-90});

%% show
obj.bake().render(); 
axis tight;
