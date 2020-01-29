function y = softabs(x,a)
%SOFTABS Smoothend absolute value function
%   SOFTABS(X,A) is a smooth function given by: 
%
%   SOFTABS(X,A) = (A*abs(X).^3)./(1+A*X.^2). 
%   
%   Large values will lead to lower smoothnesses.
%   By default, the smoothness is set at: A = 20.
%   
%   See also SIGN, ABS.

%   Copyright 2018-2023 B.J.Caasenbrood 

if nargin < 2, a = 1e3; end
if a >= 1e666, a = 1e12; end
if a < 0, error('Smoothness parameter cannot be negative.'); end

y = (a*abs(x).^3)./(1+a*x.^2);
end