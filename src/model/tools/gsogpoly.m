function [Q, R] = gsogpoly(X)
% Gram-Schmidt orthogonalization
% Written by Mo Chen (sth4nth@gmail.com).
[d,n] = size(X);
x = linspace(0,1,d).';
m = min(d,n);
R = eye(m,n);
Q = zeros(d,m);

for ii = 1:m
    R = 0;

    for jj = 1:(ii-1)
        D = (trapz(x,Q(:,jj).*Q(:,jj)));
        R = R + (trapz(x,X(:,ii).*Q(:,jj))/D)*Q(:,jj);
    end
    
    Q(:,ii) = X(:,ii) - R;
    Q(:,ii) = Q(:,ii)/sqrt(trapz(x,Q(:,ii).*Q(:,ii)));
end
