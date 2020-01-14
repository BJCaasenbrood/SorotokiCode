function d = dEllipse(P,xc,yc,a,b)
d = sqrt((b^2)*(P(:,1)-xc).^2+((a^2)*P(:,2)-yc).^2)-sqrt((a^2)*(b^2));
d=[d,d];