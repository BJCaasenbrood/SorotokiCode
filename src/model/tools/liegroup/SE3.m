function Y = SE3(X,varargin)

singleFrame = true;
if numel(X) == 7
    q = X(1:4,1);     % quaterions
    p = X(5:7,1);     % position
    R = quat2rot(q);  % Rot-mat SO(3)
elseif numel(X) == 4
    R = quat2rot(X);  % Rot-mat SO(3)    
    p = varargin{1};
elseif numel(X) == 9
    R = X;
    if ~isempty(varargin)
       p = varargin{1};
    else
       p = zeros(3,1); 
    end
else
    singleFrame = false;
end

if singleFrame
    Y          = zeros(4);
    Y(4,4)     = 1;
    Y(1:3,1:3) = R;
    Y(1:3,4)   = p;
else
    Y = zeros(4,4,size(p,1));
    Y(4,4,:)     = 1;
    Y(1:3,1:3,:) = R;
    Y(1:3,4,:)   = p;
end

end

