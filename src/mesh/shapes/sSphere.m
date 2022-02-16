function sdf = sSphere(xc,yc,zc,r)
if nargin < 2
    r = xc; 
    xc = 0;
    yc = 0;
    zc = 0;
end

sdf = Sdf(@(P) sdfSphere(P,xc,yc,zc,r));

% % r = 2*r;

sdf.BdBox = [xc-r-1e-6,xc+r+1e-6,....
             yc-r-1e-6,yc+r+1e-6,...
             zc-r-1e-6,zc+r+1e-6];
         
[sdf.Node,sdf.Element] = generateNodeSet(xc,yc,zc,r,30);
         
end

function d = sdfSphere(P,xc,yc,zc,r)
d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2+(P(:,3)-zc).^2)-r;
d = [d,d];
end

function [V,F] = generateNodeSet(xc,yc,zc,r,N)
th = linspace(-pi,pi,N).';

x1 = r*cos(th) + xc;
y1 = r*sin(th) + yc;
z1 = x1*0;

x2 = r*sin(th) + xc;
y2 = x1*0;
z2 = r*cos(th) + zc;

x3 = x1*0;
y3 = r*sin(th) + xc;
z3 = r*cos(th) + zc;

V = [x1,y1,z1;x2,y2,z2;x3,y3,z3];
S = [(1:N-1).',(2:N).'];
F = [S;S+N;S+2*N];
%F = [1:N;N+1:2*N;2*N+1:3*N];
end