function Y = admap(X,varargin)

if size(X,1) == 7
    q = R7(1:4,1);     % quaterions
    p = R7(5:7,1);     % position
    R = quat2rot(q);   % Rot-mat SO(3)
    
    S = skew(p);
    Y = zeros(6);
    Y(1:3,1:3) = R;
    Y(4:6,4:6) = R;
    Y(4:6,1:3) = S*R;
    
elseif numel(X) == 9
    R = X;
    if ~isempty(varargin)
       p = varargin{1};
    else
       p = zeros(3,1); 
    end
    
    S = skew(p);
    Y = zeros(6);
    Y(1:3,1:3) = R;
    Y(4:6,4:6) = R;
    Y(4:6,1:3) = S*R;
    
elseif numel(X) == 6
    W = X(1:3);
    U = X(4:6); 
    Wh = skew(W); 
    Uh = skew(U);
    Y = zeros(6);
    Y(1:3,1:3) = Wh;
    Y(4:6,4:6) = Wh;
    Y(4:6,1:3) = Uh;
    
elseif numel(X) == 16    
    R = X(1:3,1:3);
    p = X(1:3,4);
    
    S = skew(p);
    Y = zeros(6);
    Y(1:3,1:3) = R;
    Y(4:6,4:6) = R;
    Y(4:6,1:3) = S*R;

end

end

%--------------------------------------------------------------------------
function y = skew(v)
y = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
end