clr;
%% loading .stl file
sdf = SDF;

%% dual contouring
[F,V] = DualContour(SDF,30);

obj = Gmodel(V,F);
%obj = Gmodel(sdf,'Quality',120);
obj.bake.render

function y = SDF
    C1 = sCube(0,1,0,1,0,1);
    S1 = sSphere(0,0,1,.5);
    S2 = sSphere(0,0,0.5,1);
    
    y = (C1 - S1)/S2;
end