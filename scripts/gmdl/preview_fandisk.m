clr;
%% model
obj = Gmodel('Fandisk.stl');

%% set texture
obj.Texture = metal;
obj.Emission = [0.15 0.75 0.15];

obj.set('SSSPower',1.3,...
        'AORadius',0.6,...
        'AOPower',5,...
        'SSS',false);

obj.bake();
    
%% set view
obj.render(); 
obj.update();
axis tight;
