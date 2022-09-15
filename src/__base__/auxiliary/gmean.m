function y = gmean(x,varargin)
%GMEAN Generalized mean.
% A generalized mean, also known as power mean, Holder mean or Kolmogorov-
% Negumo function of the mean, is an abstraction of the Pythagorean means
% included harmonic, geometric, and arithmetic mean.
% It is defined as,
%
%              Mk = [1/n(x1^k + x2^k + ... + xn^k)]^1/k
%
% where: k is indicator power for the desired mean (-1 = harmonic mean;
% 0 = geometric mean; 1 = arithmetic mean;2 = root mean square).
% Although it is not possible to put k = 0 directly but, according to the
% L辿opital痴 theorem,  the limit as k tends to zero exists,
%
%              Mk = lim k->0 [1/n(x1^k + x2^k + ... + xn^k)]^1/k
%                 = (x1x2 ... xn)^1/k
%
% Syntax: function y = gmean(x)
n  = length(x);
if isempty(varargin)
    k = 0;
else
    k = varargin{1};
end

if (k == -1) || (k == 1) || (k == 2)
    xk = x.^k;
    y = (sum(xk)/n)^(1/k);
elseif (k == 0)
    lx = log(x);
    y = exp(sum(lx)/n);
end
return,