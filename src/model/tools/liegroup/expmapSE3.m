function Y = expmapSE3(X)

if numel(X) == 6
    S = isom(X(1:3));
    U = X(4:6);
elseif numel(X) == 16
    S = X(1:3,1:3);
    U = X(1:3,4);
end


Y = zeros(4);
Y(4,4)     = 1;
Y(1:3,1:3) = expmapSO3(S);
Y(1:3,4)   = (tmapSO3(S)).'*U;

end
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
function Y = tmapSO3(X)

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
    
    Y = eye(3) - b*S + (1/(t*t))*(1-a)*S*S;
else
    Y = eye(3);
end

end
%--------------------------------------------------------------------------
function y = isom(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; 
     x3, 0, -x1; 
     -x2, x1, 0];
end
%--------------------------------------------------------------------------
function y = isom_(x)
y = [x(3,2);x(1,3);x(2,1)];
end
