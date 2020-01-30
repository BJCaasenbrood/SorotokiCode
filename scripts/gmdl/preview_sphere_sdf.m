clc; clear; close all;

%% model
Dist = @(X) dSphere(X,0,0,0,1);
BdBox = [-1,1,-1,1,-1,1];

obj = Gmodel(Dist,BdBox);

%% set texture
obj.Texture = copper;
obj.bake();

%% show
obj.render(); 