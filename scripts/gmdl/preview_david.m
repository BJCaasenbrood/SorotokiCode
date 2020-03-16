clc; clear; close all;

%% model
obj = Gmodel('Cube.stl');

%% set texture
obj.Texture = 1.25*base;
% obj.Emission = lime;
% 
% obj.set('AOPower',1.5,...
%         'AORadius',0.25);

obj = Blender(obj,'Curve',{'PCC',2,2});
obj = Blender(obj,'Twist',{'z',15});

obj = obj.bake();

%% show
obj = obj.render();
