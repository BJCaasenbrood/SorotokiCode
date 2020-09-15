function varargout = cStitch(x,N)

%y = @(x) 0.5*((1-sign(x-0.5)).*c1(2*x) + (1+sign(x-0.5)).*c2(2*x - 1));
W = linspace(0,1,N+1);
for i = 1:N
    varargout{i} = 
end

end

