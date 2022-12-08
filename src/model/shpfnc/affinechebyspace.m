function [Y,X] = affinechebyspace(N,M,C,varargin)
% loop over functional space
X0 = linspace(0,1,N).';
Id = floor(linspace(1,N,C+1)); Id(1) = 0;
Nd = diff(Id);
Y  = zeros(N,C*M);

for kk = 1:C
for ii = 1:M
   W  = Nd(kk);
   X  = linspace(0,1,W).';
   Y((Id(kk)) + (1:W),(kk-1)*M + ii) = chebyshev(X,ii-1); % chebyshev
end
end

% ensure its orthonormal (gram–schmidt)
if isempty(varargin)
   Y = gsogpoly(Y,X0);
end
end

