clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture and render
obj.set('Texture',grey,'AO',true,'AOPower',5.0,'AORadius',0.2);
obj.bake().render();

view(0,15); obj.update();

%% set AO map object
obj_ = obj.copy('Translate',{'x',120});
obj_.render().showMap('AO');
view(0,15); axis tight;
