clc; clear; close all;

%% model
Dist = @(X) dSphere(X,0,0,0,1);
BdBox = [-1,1,-1,1,-1,1];

%obj = Gmodel(Dist,BdBox);
obj = Gmodel('Manifold.stl');

%obj = Blender(obj,'Scale',{'z',1.5});
%obj = Blender(obj,'Rotate',{'3D',10,20,0});

%% set texture
%obj.set('TextureStretch',.75);
obj.Texture = 1.25*grey;
obj.bake();

%% show
obj.render(); 