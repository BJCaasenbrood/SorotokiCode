clc; clear; close all;

%% model
obj = Gmodel('Astroid.stl');

%% set texture
obj.Texture = green;
% obj.Emission = lime;
% 
% obj.set('AOPower',1.5,...
%         'AORadius',0.25);

obj = obj.bake();

%% show
obj = obj.render();
