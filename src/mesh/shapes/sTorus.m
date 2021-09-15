function sdf = sTorus(xc,yc,zc,rc,rt)
if nargin < 2
    r = xc; 
    xc = 0;
    yc = 0;
    zc = 0;
end

sdf = Sdf(@(P) sdfTorus(P,xc,yc,zc,rc,rt));

sdf.BdBox = [xc-rc-rt-1e-6,xc+rc+rt+1e-6,....
             yc-rc-rt-1e-6,yc+rc+rt+1e-6,...
             zc-rc-rt-1e-6,zc+rc+rt+1e-6];
end

function d = sdfTorus(P,xc,yc,zc,rc,rt)
d = sqrt((rc-sqrt((P(:,1)-yc).^2+(P(:,2)-zc).^2)).^2 + (P(:,3)-xc).^2) - rt;
d=[d,d];
end