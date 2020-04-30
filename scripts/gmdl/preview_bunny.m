clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.Texture = greys;
obj.Emission = [0.75 0.75 0.75];

obj.set('SSSPower',0.50,...
        'AORadius',0.65,...
        'SSS',true);
    
obj.bake();

%% set view
obj.render(); 
obj.ground();

obj.showAO();
