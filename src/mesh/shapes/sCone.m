function sdf = sCone(xc,yc,zc,r,h)
if nargin < 2
    r = xc; 
    h = xc;
    xc = 0;
    yc = 0;
    zc = 0;
end

sdf = Sdf(@(P) sdfCone(P,xc,yc,zc,r,h));

sdf.BdBox = [xc-r-1e-6,xc+r+1e-6,....
             yc-r-1e-6,yc+r+1e-6,...
             zc-1e-6,zc+h+1e-6];
end

function d = sdfCone(P,xc,yc,zc,R,H)
C = R/H;
d1 = (((P(:,1)-xc).^2+(P(:,2)-yc).^2)/C^2 - (P(:,3)-zc).^2);
d3 = -(P(:,3))-H;

d = max([max([d1,d3],[],2),P(:,3)],[],2);
d = [d,d];
end
