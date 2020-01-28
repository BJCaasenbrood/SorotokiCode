function y = append(vec,x)

if size(vec,2) == 1, y = [vec;x];
else, y = [vec,x]; end

end

