function [Y,X] = pwlspace(N,M,varargin)
% loop over functional space
X  = linspace(0,1,N).';
Y0 = zeros(N,M);

for ii = 1:M
   Y0(:,ii) = pwl(X,ii,M); %chebyshev(X/L,ii-1); 
end

% ensure its orthonormal (gram–schmidt)
if isempty(varargin)
    Y = gsogpoly(Y0,X);
else
    Y = Y0;
end
end

%%
function y = pwl(x,n,c)
% n = index
% c = number of cuts
if n > c || n < 1 
   error('Index number cannot be large that number of piece-wise segments') 
end
d = (n-1)/(c-1);    
s = (-abs((x-d)*(c-1))+1);

% s1 = (x >= (n-1)/c) & (x < (n)/c);
% s2 = (x == (n)/c);
% S = s1 | s2;

if c == 1
    y = double(x)*0 + 1;
else
    y = max(double(s),0);
end
end

