clc; clear; close all;

%% model
obj = Gmodel('SoftRoboticArm.stl');

%% set texture
obj.Texture = grey;
obj = obj.bake();

%% show
obj = obj.render();
view(180,-5);
