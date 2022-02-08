function Y = SE3inv(X,varargin)

if size(X,1) == 7
    q = R7(1:4,1);     % quaterions
    p = R7(5:7,1);     % position
    R = quat2rot(q);   % Rot-mat SO(3)
elseif numel(X) == 9
    R = X;
    if ~isempty(varargin)
       p = varargin{1};
       p = p(:);
    else
       p = zeros(3,1); 
    end
end

Y          = zeros(4);
Y(4,4)     = 1;
Y(1:3,1:3) = R.';
Y(1:3,4)   = -R.'*p;

end

