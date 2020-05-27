clr;
%% generate mesh from sdf
t = 0.6;    % thickness
f = 5;      % frequency

Dist = @(X) SDF(X,t,f);
obj = Gmodel(Dist,domain(-1.1,1.1,3),'Quality',85);

%% set texture
obj.Texture = grey;
    
%% bake and render
obj.bake().render(); 

function Dist = SDF(X,t,f)
%G1 = dGyroid(X,t,f);
G1 = dSchwarzP(X,t,f);
C1 = dSphere(X,0,0,0,1.1);
Dist = dIntersect(C1,G1);
end