clc; clear; close all;

%% model
Dist = @(X) dTorus(X,0,0,0,2,1);
BdBox = [-3,3,-3,3,-3,3];

obj = Gmodel(Dist,BdBox);

%% set texture
obj.Texture = metal;
obj.bake();

%% show
obj.render(); 