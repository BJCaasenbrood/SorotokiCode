clr;
%% preview
obj = Gmodel('Bunny.stl');

%% set texture
obj.set('Texture',grey,'Emission',[0.70 0.70 0.70],...
        'SSS',true,'SSSPower',1.40,'SSSRadius',0.2);
        
obj.bake().render();

view(0,15); 
obj.update();

%% set AO map object
obj_ = obj.copy('Translate',{'x',120});
obj_.render().showMap('SSS');
view(0,15); axis tight;