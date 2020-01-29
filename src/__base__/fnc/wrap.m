function y = wrap(x,y,map)

x = mod(x,map(2)-map(1));
xlow = min(x); 
xupp = max(x);

end

