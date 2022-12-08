function sdf = sCylinder(xc,yc,z1,z2,r)
if nargin < 2
    r  = xc;
    z2 = 2*xc;
    xc = 0;
    yc = 0;
    z1 = 0;
end

sdf = Sdf(@(P) sdfCylinder(P,xc,yc,z1,z2,r));

sdf.BdBox = [xc-r-1e-6,xc+r+1e-6,....
             yc-r-1e-6,yc+r+1e-6,...
             z1-1e-6,z2+1e-6];
         
sdf.Center = sdf.centerofmass();         
         
[sdf.Node,sdf.Element] = generateNodeSet(xc,yc,z1,z2,r,30);
         
end

function d = sdfCylinder(P,xc,yc,z1,z2,r)
d1 = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
d2 = z1-P(:,3);
d3 = P(:,3)-z2;

F1 = sqrt(d1.^2 + d2.^2);
F2 = sqrt(d1.^2 + d3.^2);
F3 = max([d1,d2,d3],[],2);
L1 = (d1 > 0 & d2 >0);
L2 = (d1 > 0 & d3 >0);
L3 = ~(L1 + L2);

d = F1.*L1 + F2.*L2 + F3.*L3;

end

%-------------------------------------------------------------------------%
function [V,F] = generateNodeSet(xc,yc,zc1,zc2,r,N)
th = linspace(-pi,pi,N).';
tt = linspace(-pi,pi-pi/2,4).';
zh = linspace(zc1,zc2,N).';
nn = 0*th;

x1 = r*cos(th) + xc;
y1 = r*sin(th) + yc;
z1 = zc1 + nn;

x2 = r*cos(th) + xc;
y2 = r*sin(th) + yc;
z2 = zc2 + nn;

x3 = r*cos(tt(1)) + xc + nn;
y3 = r*sin(tt(1)) + yc + nn;
z3 = zh;

x4 = r*cos(tt(2)) + xc + nn;
y4 = r*sin(tt(2)) + yc + nn;
z4 = zh;

x5 = r*cos(tt(3)) + xc + nn;
y5 = r*sin(tt(3)) + yc + nn;
z5 = zh;

x6 = r*cos(tt(4)) + xc + nn;
y6 = r*sin(tt(4)) + yc + nn;
z6 = zh;

V = [x1,y1,z1;x2,y2,z2;x3,y3,z3;...
     x4,y4,z4;x5,y5,z5;x6,y6,z6];
 
S = [(1:N-1).',(2:N).'];
F = [S;S+N;S+2*N;S+3*N;S+4*N;S+5*N];
end