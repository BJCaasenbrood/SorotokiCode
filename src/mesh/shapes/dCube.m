function d = dCube(P,x1,x2,y1,y2,z1,z2)
d = [x1-P(:,1) + 1e-12, P(:,1)-x2, y1-P(:,2)+1e-12, P(:,2)-y2, z1-P(:,3)+1e-12, P(:,3)-z2];
d = [d,max(d,[],2)];
