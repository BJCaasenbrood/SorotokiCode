function d = dTorus(P,xc,yc,zc,c,r)
d = sqrt((c-sqrt((P(:,2)-yc).^2+(P(:,3)-zc).^2)).^2 + (P(:,1)-xc).^2) - r;
d=[d,d];
