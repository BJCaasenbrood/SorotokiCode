clc; clear; close all;

%% model
sdf = sTorus(0,0,0,2,1);
obj = Gmodel(sdf);

%% set texture
obj.set('Texture',softmath);
obj.bake().render(); 
   
