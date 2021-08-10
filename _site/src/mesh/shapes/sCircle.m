function sdf = sCircle(xc,yc,r)
if nargin < 2
    r = xc; 
    xc = 0;
    yc = 0;
end

sdf = SDF(@(P) dCircle(P,xc,yc,r));
sdf.BdBox = [xc-r-1e-6,xc+r+1e-6,yc-r-1e-6,yc+r+1e-6];
end

function d = dCircle(P,xc,yc,r)
d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
d=[d,d];
end