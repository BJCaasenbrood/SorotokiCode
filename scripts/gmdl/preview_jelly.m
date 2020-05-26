clc; clear; close all;

%% model
obj = Gmodel('David.stl');

%% set texture
obj.Texture = greenglass;
obj.Emission = lime;

obj.set('SSSPower',.75,'SSSRadius',0.3,'SSS',true);

%% set visual
obj.bake().render().update(); 