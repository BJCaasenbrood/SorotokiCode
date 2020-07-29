clc; clear; close all;

%% model
obj = Gmodel('SoftRoboticArm.stl');

%% set texture
obj.Texture = studioclay;
obj.Emission = [1.0,0.95,0.75];
obj.set('SSSPower',0.25,'SSSRadius',0.3,'SSS',true);
obj.set('AOPower',0.95,'AORadius',0.3,'AO',true);

%% set visual
obj.bake().render().update(); 
