function sdf = sLine(x1,x2,y1,y2)
if nargin < 2
    if nume(x1) == 4
       y2 = x1(4);
       y1 = x1(3);
       x2 = x1(2);
       x1 = x1(1);
    end
end

sdf = SDF(@(P) sdfLine(P,x1,x2,y1,y2));
sdf.BdBox = [x1,x2,y1,y2];
end

function d = sdfLine(P,x1,x2,y1,y2)
a = [x2-x1,y2-y1]; 
a = a/norm(a);
b = [P(:,1)-x1,P(:,2)-y1];
d = b(:,1)*a(2) - b(:,2)*a(1);
d = [d,d];
end