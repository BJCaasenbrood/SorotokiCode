%------------------------------------------------- SMOOTH ABSOLUTE FUNCTION
function y = softabs(x,a)
if nargin < 2, a = 20; end
y = (a*abs(x).^3)./(1+a*x.^2);
end