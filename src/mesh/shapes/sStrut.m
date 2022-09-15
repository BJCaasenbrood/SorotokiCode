function sdf = sStrut(x1,x2,y1,y2,T)
if nargin < 2
    if nume(x1) == 4
       y2 = x1(4);
       y1 = x1(3);
       x2 = x1(2);
       x1 = x1(1);
    end
end

sdf = Sdf(@(P) sdfLine(P,x1,x2,y1,y2,T));
sdf.BdBox = [x1-T,x2+T,y1-T,y2+T];
end

function d = sdfLine(P,x1,x2,y1,y2,T)
a = [x2-x1,y2-y1]; 
a = a/norm(a);
b = [P(:,1)-x1,P(:,2)-y1];
d1 = abs(b(:,1)*a(2) - b(:,2)*a(1))-T;

th = atan2(y2-y1,x2-x1);
r = [cos(-pi/2+th), sin(-pi/2+th);
     cos(+pi/2+th), sin(+pi/2+th)];

X11 = x1;
X21 = x1 + r(1,1);
Y11 = y1;
Y21 = y1 + r(1,2);

X12 = x2;
X22 = x2 + r(2,1);
Y12 = y2;
Y22 = y2 + r(2,2);

a = [X21-X11,Y21-Y11]; 
a = a/norm(a);
b = [P(:,1)-X11,P(:,2)-Y11];
d2 = b(:,1)*a(2) - b(:,2)*a(1);

d1=max(d1(:,end),d2(:,end));

a = [X22-X12,Y22-Y12]; 
a = a/norm(a);
b = [P(:,1)-X12,P(:,2)-Y12];
d2 = b(:,1)*a(2) - b(:,2)*a(1);

d=max(d1(:,end),d2(:,end));

c1 = sqrt((P(:,1)-x1).^2+(P(:,2)-y1).^2)-T;
c2 = sqrt((P(:,1)-x2).^2+(P(:,2)-y2).^2)-T;

d = min(min(d(:,end),c2(:,end)),c1(:,end));

d = [d,d];
end