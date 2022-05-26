clc; clear; close all;

%% model
sdf = sSphere(1);
obj = Gmodel(sdf,'Quality',10);

%% set texture
obj.Texture = softmath;
% or
obj.set('Texture',softmath);
obj.bake().render(); 
   