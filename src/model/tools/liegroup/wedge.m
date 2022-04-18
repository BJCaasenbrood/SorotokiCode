function Y = wedge(X)

if numel(X) == 16
    S = X(1:3,1:3);
    U = X(1:3,4);
    
    Y = zeros(6,1);
    Y(1:3) = isom_(S);
    Y(4:6) = U;
elseif numel(X) == 9
    Y = isom_(X);
end

end

%--------------------------------------------------------------------------
function y = isom_(x)
y = [x(3,2);x(1,3);x(2,1)];
end

