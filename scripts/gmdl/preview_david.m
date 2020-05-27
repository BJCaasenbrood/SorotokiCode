clc; clear; close all;

%% model
obj = Gmodel('David.stl');

%% set texture
obj.Texture = porcelain;
%obj.set('AO',true,'AOPower',3.0,'AORadius',0.25);

obj = obj.bake().render();
