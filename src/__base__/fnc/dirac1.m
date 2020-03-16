function Y = dirac1(n,X)
%DIRAC  Delta function.
%    DIRAC(X) is zero for all X, except X == 0 where it is infinite.
%    DIRAC(n,X) is the nth derivative of DIRAC(X).
%    DIRAC(X) is not a function in the strict sense, but rather a
%    distribution with int(dirac(x-a)*f(x),-inf,inf) = f(a) and
%    diff(heaviside(x),x) = dirac(x).
%    See also HEAVISIDE.

%   Copyright 1993-2014 The MathWorks, Inc.

if nargin == 1
    X = n;
    n = 0;
end
if ~isequal(round(n), n) || any(~isfinite(n(:))) || ~isreal(n) || any(n(:) < 0)
  error(message('symbolic:dirac:ExpectingNonNegativeInteger1'));
end
if ~isequal(size(X),size(n)) && ~isscalar(X) && ~isscalar(n)
  error(message('symbolic:sym:ScalarOrSameSize'));
end
% we want to distinguish 0 and -0 in X, scalar expansion throws that away
Xinv = -1./X;
% scalar expansion, make n and X same size - simplifies code for Y below
n = n+zeros(size(X));
X = X+zeros(size(n));
Xinv = Xinv+zeros(size(X));
Y = zeros(size(X),'like',X);

% dirac(n, +0) = (-1)^n*Inf, dirac(n, -0) = Inf
Y(X == 0) = Xinv(X==0);
Y(X == 0 & round(n/2) == n/2) = 1;

Y(imag(X) ~= 0) = 0;
Y(isnan(X)) = 0;
end
