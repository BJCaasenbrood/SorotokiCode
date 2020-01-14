function d = dSphere(P,xc,yc,zc,r)
d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2+(P(:,3)-zc).^2)-r;
d = [d,d];
