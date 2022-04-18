function Y = hat(X)

if numel(X) == 6
    S = skew(X(1:3));
    U = X(4:6);
    
    Y = zeros(4,4);
    Y(1:3,1:3) = S;
    Y(1:3,4)   = U;
    
elseif numel(X) == 3
    Y = skew(X);
end
    
end

%--------------------------------------------------------------------------
function y = skew(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; 
     x3, 0, -x1; 
     -x2, x1, 0];
end

