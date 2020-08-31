clc; clear; close all;

%% model
Dist = @(x) SDF(x);
obj = Gmodel(Dist,domain(0,1,3),'Quality',150);

%% set texture
obj.Texture = clean;

%% show
obj.bake().render(); 

function Dist = SDF(x)
Dist = dCube(x,0,1,0,1,0,0.01);
end