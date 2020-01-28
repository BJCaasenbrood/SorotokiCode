clc; clear; close all;

%% model
obj = Gmodel('SoftRoboticArm.stl');

%% set texture
obj.Texture = ImageColormix(retro,[49,200,255]);
%obj.Emission = lime;

% obj.set('AOPower',1.5,...
%         'AORadius',0.25);

obj = obj.bake();

%% show
obj = obj.render();
