function Y = admapinv(X,varargin)

if size(X,1) == 7
    q = R7(1:4,1);     % quaterions
    p = R7(5:7,1);     % position
    R = quat2rot(q);   % Rot-mat SO(3)
    
elseif numel(X) == 9
    R = X;
    if ~isempty(varargin)
       p = varargin{1};
    else
       p = zeros(3,1); 
    end
elseif numel(X) == 16
    R = X(1:3,1:3);
    p = X(1:3,4);    
end

Rt = R.';
S  = skew(p);
Y = zeros(6);
Y(1:3,1:3) = Rt;
Y(4:6,4:6) = Rt;
Y(4:6,1:3) = Rt*S.';

end

function y = skew(x)
    x1 = x(1); x2 = x(2); x3 = x(3);
    y = [0, -x3, x2; x3, 0, -x1; -x2, x1, 0];
end
