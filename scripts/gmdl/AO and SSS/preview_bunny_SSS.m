clr;
%% preview
obj = Gmodel('Bunny.stl');

%% set texture
obj.set('Texture',grey,'Emission',[0.70 0.70 0.70],...
        'SSS',true,'SSSPower',1.70,'SSSRadius',0.2);
        
obj.bake().render().update();

%% set AO map object
obj_ = obj.copy('Translate',{'y',100});
obj_.render().showMap('SSS');
view(90,15); axis('tight');