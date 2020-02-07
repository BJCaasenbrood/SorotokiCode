clr;
%% set signed distance function
msh = Mesh(@(x) SDF(x,0),'BdBox',[-1 3 -1 3]);
figure(101), subplot(2,2,1); zoom(1.33);
msh.showSDF;

msh = Mesh(@(x) SDF(x,1),'BdBox',[-2 2 -1 3]);
figure(101), subplot(2,2,2); zoom(1.33);
msh.showSDF;

msh = Mesh(@(x) SDF(x,2),'BdBox',[-1.5 2.5 -1 3]);
figure(101), subplot(2,2,3); zoom(1.33);
msh.showSDF;

msh = Mesh(@(x) SDF(x,3),'BdBox',[-1 3 -1 3]);
figure(101), subplot(2,2,4); zoom(1.33);
msh.showSDF;

function Dist = SDF(x,b)
R = dRectangle(x,0,2,0,2);
C = dCircle(x,0,1,1.0);
if b == 0, Dist = R;
elseif b == 1, Dist = C;
elseif b == 2, Dist = dUnion(R,C);
elseif b == 3, Dist = dDiff(R,C);
end
end

