function [L,dl] = arclength(X)
%[~,dX] = gradient(X);
dX = diff(X);
dl = sqrt(sum(dX.^2,2));
L  = sum(dl);
end

