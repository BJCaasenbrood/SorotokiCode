function d = dCuboid(P,a,b,c)
d = [max([P(:,1).^2 - a.^2, P(:,2).^2 - b.^2],[],2),P(:,3).^2 - c^2];
d = [d,max(d,[],2)];
end