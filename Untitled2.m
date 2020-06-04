clr;
%% generating graphical models
obj0 = Gmodel('Bunny.stl');
obj1 = Gmodel(@(x) SDF(x),domain(-0.1,1.1,3));

%% rendering models
figure(101);
subplot(1,2,1); obj0 = obj0.bake().render(); 
subplot(1,2,2); obj1 = obj1.bake().render(); 

%% signed distance fucntion (3D)
function Dist = SDF(x)
    C1 = dCube(x,0,1,0,1,0,1);
    S1 = dSphere(x,0,0,1,.5);
    S2 = dSphere(x,0,0,0.5,1);
    Dist = dIntersect(dDiff(C1,S1),S2);
end