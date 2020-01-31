clc; clear; close all;

%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.Texture = base;
obj.Emission = [0.65 0.65 0.65];

obj.set('SSSPower',1.2,...
        'AORadius',0.3,...
        'AOPower',1.0,...
        'SubSurfaceScattering',true);

obj.bake();
    
%% set view
obj.render(); 
view(260,20); 
obj.update();
