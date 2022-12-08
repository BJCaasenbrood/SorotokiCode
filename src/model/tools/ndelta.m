function d = ndelta(x,h)
% nascent delta function
d = (1/(sqrt(pi)))*exp(-(x.^2)/h^2);
end
