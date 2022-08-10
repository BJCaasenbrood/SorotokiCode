function y = srsin(x,e)
y = ((sqrt(1 + e.^2)))*sin(x)./(sqrt(sin(x).^2 + e.^2));
end

