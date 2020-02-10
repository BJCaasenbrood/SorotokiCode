clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.Texture = grey;
obj.Emission = [0.75 0.75 0.75];

obj.set('SSSPower',1.0,...
        'AORadius',0.2,...
        'AOPower',5,...
        'AmbientOcclusion',true);

obj.bake();
    
%% set view
obj.render(); 
view(240,20); 
obj.update();
