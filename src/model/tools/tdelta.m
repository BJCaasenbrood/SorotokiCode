function y = tdelta(x,h)
y = max(1-abs((2*x)/(h) - (0.5)/(h)),0);
end

