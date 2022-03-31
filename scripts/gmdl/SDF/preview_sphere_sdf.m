clc; clear; close all;

%% model
sdf = sSphere(0,0,0,1);
obj = Gmodel(sdf,'Quality',150);

%% set texture
obj.Texture = clean;

%% show
obj.bake().render(); 