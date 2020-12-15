function sdf = sRectangle(x1,x2,y1,y2)
if nargin < 2
    a = x1;
    x1 = -a;
    x2 = a;
    y1 = -a;
    y2 = a;
end

eps = 5e-2*norm([x1,x2,y1,y2]);

sdf = SDF(@(P) dRectangle(P,x1,x2,y1,y2));
sdf.BdBox = [x1-eps,x2+eps,y1-eps,y2+eps];
end

function d = dRectangle(P,x1,x2,y1,y2)
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2];
d = [d,max(d,[],2)];
end