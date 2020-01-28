clc; clear; close all;

%% model
obj = Gmodel('David.stl');

%% set texture
obj.Texture = highway;
obj.Emission = lime;

obj.set('AOPower',1.5,...
        'AORadius',0.25);

obj = obj.bake();

%% show
obj = obj.render();
