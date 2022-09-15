function X = carc(p1,p2,r)
p = [p1(1),p2(1);p1(2),p2(2)]; % p - this is our points p1 & p2, p =[x1,x2;y1,y2];
%r = 7;          % r - radius
% Finding the coordinates of the centers of circles
% xy = [x1, y1; x2, y2]
a = sym('a',[2,1],'real');
eqs = [1,1]*(p - repmat(a(:),1,2)).^2 - abs(r)^2;
sol = vpasolve(eqs,a);
ss = struct2cell(sol);
xy = double([ss{:}]);
% example: The arc of a circle with center at xy(1,:)
v = xy(1,:);
p1 = p - v(:);
alp = atand(p1(2,:)./p1(1,:));
alp = alp + 180*(p1(1,:) < 0 & p1(2,:) > 0) - 180*(p1(1,:) < 0 & p1(2,:) < 0);
asort = sort(alp);
if sign(r) == 1
    phi = linspace(asort(1),asort(2),100)';
else
    dd = asort(2) - asort(1);
    phi = linspace(-asort(1),asort(2)-(360-dd),100)';
end
X = [ abs(r)*cosd(phi(:)) + v(1), abs(r)*sind(phi(:)) + v(2)];

end


