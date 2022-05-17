clc; clear; close all;

%% model
sdf = sSphere(0,0,0,1);
obj = Gmodel(sdf,'Quality',20);

%% set texture
obj.set('Texture',softmath);
obj.bake().render(); 
   