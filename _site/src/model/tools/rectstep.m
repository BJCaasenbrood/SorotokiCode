function y = rectstep(x,a,b)
y = zeros(numel(x),1);
y(x>=a & x<=b) = 1;
end

