clc; clear; close all;

%% model
obj = Gmodel('Dragon.stl');

%% set texture
obj.Texture = wax;
obj.Emission = [0.25 0.85 0.55];

obj.set('SSSPower',0.72,...
        'AORadius',0.3,...
        'AOPower',0.40,...
        'SubSurfaceScattering',true);

obj.bake();
    
%% set view
obj.render(); 
view(260,20); 
obj.update();
