clc; clear; close all;

%% model
Dist = @(X) dTorus(X,0,0,0,2,1);
BdBox = [-3,3,-3,3,-3,3];

obj = Gmodel(Dist,BdBox);

%% set texture
obj.Texture = jade;
obj.Emission = [0.15 0.75 0.15];

obj.set('SSSPower',1.3,'SSS',true);
obj.bake();
    

%% show
obj.render(); 