function d = dndelta(x,h)
% nascent delta function

% if abs(x) >= h, d = 0;
% elseif (x >= -h && x <= 0), d = (x/h^2) + 1/h; 
% else, d = -x/h^2 + 1/h;
% end

d = -(2*x*exp(-x^2/h^2))/(h^2*pi^(1/2));

end
