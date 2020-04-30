clr;
%% generate mesh from sdf
t = 0.5;    % thickness
f = 4;      % frequency

Dist = @(X) SDF(X,t,f);
BdBox = domain(-1.1,1.1,3);
obj = Gmodel(Dist,BdBox,'Quality',50);

%% set texture
obj.Texture = grey;
obj.bake();

%% show
obj.render(); 

function Dist = SDF(X,t,f)
%G1 = dGyroid(X,t,f);
G1 = dSchwarzP(X,t,f);
R1 = dCube(X,-1,1,-1,1,-1,1);
Dist = dIntersect(R1,G1);
end