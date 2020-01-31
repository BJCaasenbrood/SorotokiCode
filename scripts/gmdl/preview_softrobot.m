clc; clear; close all;

%% model
obj = Gmodel('SoftRoboticArm.stl');

%% set texture
obj.Texture = base;
obj = obj.bake();

%% show
obj = obj.render();
