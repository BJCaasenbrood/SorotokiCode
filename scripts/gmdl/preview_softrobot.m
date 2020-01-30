clc; clear; close all;

%% model
obj = Gmodel('SoftRoboticArm.stl');

%% set texture
obj.Texture = green;
obj = obj.bake();

%% show
obj = obj.render();
