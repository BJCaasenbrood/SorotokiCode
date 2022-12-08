function d = dndelta(x,h)
% nascent delta function
d = -(2*x.*exp(-x.^2/h^2))/(h^2*pi^(1/2));
end
