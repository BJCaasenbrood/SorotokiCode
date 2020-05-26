clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture and render
obj.set('Texture',grey,'AO',true,'AOPower',5.0,'AORadius',0.2);
obj.bake().render().update();

%% set AO map object
obj_ = obj.copy('Translate',{'y',100});
obj_.render().showMap('AO');
view(90,15); axis tight;
