clc; clear; close all;

%% model
Dist = @(X) dSphere(X,0,0,0,1);
obj = Gmodel(Dist,domain(-1,1,3),'Quality',150);

%% set texture
obj.Texture = clean;

%% show
obj.bake().render(); 