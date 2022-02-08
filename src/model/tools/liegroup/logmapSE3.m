function Y = logmapSE3(X)

if numel(X) == 7
    q = X(1:4,1);      % quaterions
    p = X(5:7,1);      % position
    R = quat2rot(q);   % Rot-mat SO(3)
elseif numel(X) == 16
    R = X(1:3,1:3);
    p = X(1:3,4);
end

S = logmapSO3(R);
Y = zeros(4);
Y(1:3,4)   = (tmapinvSO3(S)).'*p;
Y(1:3,1:3) = S;

end
%--------------------------------------------------------------------------
function Y = tmapinvSO3(X)

if numel(X) == 3
    S = isom(X);
elseif numel(X) == 9
    S = X;
    X = isom_(S);
end

t = norm(X);
if abs(t) >= 1e-6
    a = sin(t)/t;
    b = 2*(1-cos(t))/(t*t);
    
    Y = eye(3) - 0.5*S + (1/(t*t))*(1-a/b)*S*S;
else
    Y = eye(3);
end

end
%--------------------------------------------------------------------------
function Y = logmapSO3(X)

if numel(X) == 3
    S = isom(X);
elseif numel(X) == 9
    S = X;
    X = isom_(S);
end

theta = acos(0.5*(trace(S) - 1));
alpha = 0.5*theta/sin(theta);

if ~isnan(alpha)
   Y = alpha*(S - S.'); 
else
   Y = zeros(3); 
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
