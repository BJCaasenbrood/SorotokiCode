function d = ndelta(x,h)
% nascent delta function

% if abs(x) >= h, d = 0;
% elseif (x >= -h && x <= 0), d = (x/h^2) + 1/h; 
% else, d = -x/h^2 + 1/h;
% end

d = (1/(sqrt(pi)))*exp(-(x.^2)/h^2);

end
