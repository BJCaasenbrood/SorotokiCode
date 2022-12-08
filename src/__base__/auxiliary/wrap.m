function y = wrap(x,map)
y = mod(x,map(2)-map(1))+map(1);
end

