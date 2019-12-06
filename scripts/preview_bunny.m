clc; clear; close all;

%% load model
obj = Gmodel('Bunny.stl');

%% set texture
obj.Texture = jade;
obj = obj.bake();

%% show
obj.show();