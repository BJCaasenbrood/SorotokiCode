clr;
%% model
obj = Gmodel('Soft_Module_TPU.stl');
obj.set('TextureStretch',0.95);

%% deform
%Blender(obj,'Translate',{'y',-0.4});
Blender(obj,'Curve',{'PCC+',60,0,1.2});
%Blender(obj,'Rotate',{'y',-90});

%% show
obj.bake().render(); 
axis tight;
