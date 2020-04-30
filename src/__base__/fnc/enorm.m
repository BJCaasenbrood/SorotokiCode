function [sigma,e] = enorm(y1,y2)
%ENORM Summary of this function goes here
x1 = linspace(0,1,length(y1));
x2 = linspace(0,1,length(y2));

vq = interp1(x2,y2,x1);
e = y1(:)-vq(:);
sigma = norm(e);
end

