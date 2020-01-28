clc; clear; close all;

%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.Texture = mateplastic;
obj.Emission = [0.65 0.65 0.65];

obj.set('AOPower',1.0,...
        'AORadius',0.3,...
        'SubSurfaceScattering',true);

obj.bake();

%% show
obj.render(); 

%% set view
view(260,20); obj.update();
