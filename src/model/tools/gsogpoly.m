function [Q, R] = gsogpoly(Y,varargin)
% Gram-Schmidt orthogonalization
% Written by Mo Chen (sth4nth@gmail.com).
[d,n] = size(Y);
if isempty(varargin)
    x = linspace(0,1,d).';
else
    x = varargin{1}; 
end
m = min(d,n);
R = eye(m,n);
Q = zeros(d,m);

for ii = 1:m
    R = 0;

    for jj = 1:(ii-1)
        D = (trapz(x,Q(:,jj).*Q(:,jj)));
        R = R + (trapz(x,Y(:,ii).*Q(:,jj))/D)*Q(:,jj);
    end
    
    Q(:,ii) = Y(:,ii) - R;
    Q(:,ii) = Q(:,ii)/sqrt(trapz(x,Q(:,ii).*Q(:,ii)));
end
