clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.Texture = grey;
obj.Emission = [0.75 0.75 0.75];

obj.set('SSSPower',1.0,'SSSRadius',0.3,'SSS',true,...
        'AOPower',1.5,'AORadius',0.3,'AO',true);
    
obj.bake().render(); 