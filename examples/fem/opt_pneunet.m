clr;
%% generate mesh from sdf
W = 20;   % width cell
H = 40;   % width cell
D = 1.5;  % inter distance

f = @(x) PneuNet(x,W,H,D,W);
sdf = Sdf(f,'BdBox',[0,W,0,H]);

msh = Mesh(sdf,'NElem',1200);
msh = msh.generate();

function Dist = PneuNet(P,W,H,E,T)
    R1 = dRectangle(P,0,W,0,H);
    R2 = dRectangle(P,-W/2,E,T,H+H/2);
    R3 = dRectangle(P,W-E,W+W/2,T,H+H/2);
    C1 = dCircle(P,0,T + 0.5,1);
    C2 = dCircle(P,W,T + 0.5,1);
    Dist = dDiff(dDiff(dDiff(dDiff(R1,R2),R3),C1),C2);
end