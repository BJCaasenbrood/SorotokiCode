clr;
%% model
obj = Gmodel('Pneunet.stl');
obj = obj.set('Shading','Face');

%% set texture
obj.Texture = base;
obj = obj.bake();

%% deform
obj = Blender(obj,'Rotate',{'y',90});
obj = Blender(obj,'Curve',{'PCC+',120,0,0.7});
obj = Blender(obj,'Rotate',{'y',-90});

%% show
obj = obj.render();
axis tight
