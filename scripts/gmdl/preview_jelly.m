clc; clear; close all;

%% model
obj = Gmodel('Pneunet.stl');

%% set texture
obj.Texture = redwax;
obj.Emission = 1.5*col(2);

obj.set('SSSPower',0.75,...
        'AORadius',0.3,...
        'SSS',true);

obj.bake();
    
%% set view
obj.render(); 
view(260,20); 
obj.update();
