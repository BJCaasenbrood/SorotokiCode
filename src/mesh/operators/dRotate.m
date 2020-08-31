function y = dRotate(x,rot,point)

if nargin == 3
y = x - repmat(point(:)',length(x),1);
else, y = x;
end

k = rot*(pi/180);

R = [cos(k),-sin(k);...
     sin(k), cos(k)];

y = transpose(R*y.');

if nargin == 3
y = y + repmat(point(:)',length(x),1);
end

end
%-------------------------------------------------------------------------%

