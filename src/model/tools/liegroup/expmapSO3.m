function Y = expmapSO3(X)

if numel(X) == 3
    S = isom(X);
elseif numel(X) == 9
    S = X;
    X = isom_(S);
end

t = norm(X);
if abs(t) >= 1e-6
    a = sin(t)/t;
    b = (1-cos(t))/(t*t);
    
    Y = eye(3) + a*S + b*S*S;
else
    Y = eye(3);
end

end

function y = isom(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; 
     x3, 0, -x1; 
     -x2, x1, 0];
end

function y = isom_(x)
y = [x(3,2);x(1,3);x(2,1)];
end
