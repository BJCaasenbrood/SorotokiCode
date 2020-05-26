function y = sat(x,bdn)
if nargin < 2, bdn = [-1,1]; end
y = (atan(5.*x)-atan(-5.*x))/(pi);
end

