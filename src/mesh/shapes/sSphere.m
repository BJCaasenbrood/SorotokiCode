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
end

function d = sdfSphere(P,xc,yc,zc,r)
d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2+(P(:,3)-zc).^2)-r;
d = [d,d];
end