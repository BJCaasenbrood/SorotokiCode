clr;
%% model
Dist = @(x) SDF(x);
obj = Gmodel(Dist,domain(-0.1,1.1,3),'Quality',50);

%% set texture
obj.Texture = base;

%% show
obj.bake().render(); 

function Dist = SDF(x)
C1 = dCube(x,0,1,0,1,0,1);
S1 = dSphere(x,0,0,1,.5);
S2 = dSphere(x,0,0,0.5,1);
Dist = dIntersect(dDiff(C1,S1),S2);
end