function d = dCircle(P,xc,yc,r)
d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
d=[d,d];
end