function [Y,X] = chebyspace(N,M)
% loop over functional space
X = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

