function y = bellc(x,varargin)
% BELLC produces a gaussian bell curve with color parameters sigma and mu,
%   being the standard dev. and the mean value: y = bellc(x,sigma,mu).
%   Usages:
%
%   x  = linspace(-10,10,1e3);   % sampling
%   y1 = bellc(x);               % bell curve with sigma = 1, mu = 0
%   y1 = bellc(x,5)              % bell curve with sigma = 5, mu = 0
%   y1 = bellc(x,5,-1)           % bell curve with sigma = 5, mu = -1.
%
%   Last edit: July 20, 2022.
if ~isempty(varargin)
    if numel(varargin) == 1
        sigma = varargin{1};
        mu = 0;
    elseif numel(varargin) == 2
        sigma = varargin{1};
        mu    = varargin{2};
    end
else
    sigma = 1;
    mu = 0;
end
    
beta = -1/(2*sigma^2);
y = (1/(sigma*2*pi))*exp(beta*(x-mu).^2);
end

