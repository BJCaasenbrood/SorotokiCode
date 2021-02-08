function  [a,b] = polygonArea(x,y)
%  Copyright (c) 1995 by Kirill K. Pankratov,
%	kirill@plume.mit.edu.
%	04/20/94, 05/20/95  

 % Make polygon closed ( even if it already is) .............
x = [x(:); x(1)];
y = [y(:); y(1)];

 % Calculate contour integral Int -y*dx  (same as Int x*dy).
lx = length(x);
if isfloat(x) && isfloat(y) 
a  = -(x(2:lx)-x(1:lx-1))'*(y(1:lx-1)+y(2:lx))/2;
else % mtimes can not handle ints: return double
a  = -sum((x(2:lx)-x(1:lx-1)).*(y(1:lx-1)+y(2:lx)))/2;
end

b = sign(a)==-1;

end