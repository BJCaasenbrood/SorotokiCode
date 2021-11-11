clr;
%% set signed distance function
msh = Mesh(SDF(0),'BdBox',[-1 3 -1 3]);
subplot(2,2,1); msh.showSDF;

msh = Mesh(SDF(1),'BdBox',[-2 2 -1 3]);
subplot(2,2,2); msh.showSDF;

msh = Mesh(SDF(2),'BdBox',[-1.5 2.5 -1 3]);
subplot(2,2,3); msh.showSDF; 

msh = Mesh(SDF(3),'BdBox',[-1 3 -1 3]);
subplot(2,2,4); msh.showSDF; 

function Dist = SDF(b)

R = sRectangle(0,2,0,2);
C = sCircle(0,1,1);

switch(b)
    case(0), Dist = R;
    case(1), Dist = C;
    case(2), Dist = R+C;
    case(3), Dist = R-C;
end
end

