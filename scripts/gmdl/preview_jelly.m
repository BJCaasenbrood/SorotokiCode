clc; clear; close all;

%% model
obj = Gmodel('SoftRoboticArm.stl');

%% set texture
obj.Texture = studioclay;
obj.Emission =  [1.00 0.95 0.85];

%obj.set('SSSPower',0.125,'SSSRadius',0.3,'SSS',true);
%obj.set('AOPower',0.75,'AORadius',0.3,'AO',true);

%% set visual
obj.bake().render().update(); 