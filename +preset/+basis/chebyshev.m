function Y = chebyshev(varargin)

    p = inputParser;
    addOptional(p,'n',100);
    addOptional(p,'m',3);
    parse(p,varargin{:});

    [Y, ~] = polyspace(p.Results.n, p.Results.m);
end

function [Y, X] = polyspace(N,M)
% loop over functional space
X = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(X,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
% if isempty(varargin)
%     Y = gsogpoly(Y,X);
% end
end
    
    