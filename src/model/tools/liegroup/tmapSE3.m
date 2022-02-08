% Computations can be found Sonneville, V., Cardona, A., & Brüls, O. (2014). 
% Geometrically exact beam finite element formulated on the special 
% Euclidean group SE(3). Computer Methods in Applied Mechanics and 
% Engineering, 268, 451–474. https://doi.org/10.1016/j.cma.2013.10.008
function Y = tmapSE3(X)

if numel(X) == 7
    q = X(1:4,1);      % quaterions
    p = X(5:7,1);      % position
    R = quat2rot(q);   % Rot-mat SO(3)
    
    Z = logmapSE3(SE3(R,p));
    S = isom_(Z(1:3),Z(1:3));
    U = Z(1:3,4);
    
elseif numel(X) == 16 && abs(X(4,4)-1) <= 1e-6 
    R = X(1:3,1:3);
    p = X(1:3,4);
    
    Z = logmapSE3(SE3(R,p));
    S = isom_(Z(1:3,1:3));
    U = Z(1:3,4);
    
elseif numel(X) == 16 && abs(X(4,4)) <= 1e-6 
    S = isom_(X(1:3,1:3));
    U = X(1:3,4);    
    
elseif numel(X) == 6 
    S = X(1:3);
    U = X(4:6);
end


Ts = tmapSO3(S);
Tu = touplusSE3(S,U);

Y = zeros(6);
Y(1:3,1:3) = Ts;
Y(4:6,4:6) = Ts;
Y(1:3,4:6) = Tu;

end
%--------------------------------------------------------------------------
function Y = touplusSE3(S,U)
A = isom(U); B = isom(S);
C = A*B + B*A;

t = norm(S);
if abs(t) >= 1e-6
    a = sin(t)/t;
    b = 2*(1-cos(t))/(t*t);
    
    Y = -0.5*b*A + (1/(t*t))*(1-a)*C + ...
        (1/(t*t))*(B.')*A*((b-a)*B + ...
        (0.5*b - (1/(t*t))*(3*(1-a)))*B*B);
else
    Y = -0.5*A;
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
