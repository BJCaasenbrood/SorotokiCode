function d = dCone(P,xc,yc,zc,R,H)
C = R/H;
d1 = (((P(:,1)-xc).^2+(P(:,2)-yc).^2)/C^2 - (P(:,3)-zc).^2);
d3 = -(P(:,3))-H;

d = max([max([d1,d3],[],2),P(:,3)],[],2);
