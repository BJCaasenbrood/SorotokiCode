function y = smoothstep(X)

for ii = 1:length(X)
    x = X(ii);
    if x<=0, y(ii,1) = 0;
    elseif (x>0 && x<=1), y(ii,1) = 3*x.^2 -2*x.^3;
    else, y(ii,1) = 1;
    end
end

end