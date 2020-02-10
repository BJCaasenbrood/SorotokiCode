clr;
h = tight_subplot(2,2,[0.1 0.02]);

%% set signed distance function
msh = Mesh(@(x) SDF(x,0),'BdBox',[-1 3 -1 3]);
figure(101), axes(h(1)); zoom(1.6);
msh.showSDF;

msh = Mesh(@(x) SDF(x,1),'BdBox',[-2 2 -1 3]);
figure(101), axes(h(2)); zoom(1.6);
msh.showSDF;

msh = Mesh(@(x) SDF(x,2),'BdBox',[-1.5 2.5 -1 3]);
figure(101), axes(h(3)); zoom(1.6);
msh.showSDF;

msh = Mesh(@(x) SDF(x,3),'BdBox',[-1 3 -1 3]);
figure(101), axes(h(4)); zoom(1.6);
msh.showSDF;

function Dist = SDF(x,b)
R = dRectangle(x,0,2,0,2);
C = dCircle(x,0,1,1.0);
switch(b)
    case(0), Dist = R;
    case(1), Dist = C;
    case(2), Dist = dUnion(R,C);
    case(3), Dist = dDiff(R,C);
end
end

