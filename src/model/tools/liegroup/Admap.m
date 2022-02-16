function Ad = Admap(X,varargin)

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
end


S = isom(p);
Ad = zeros(6);
Ad(1:3,1:3) = R;
Ad(4:6,4:6) = R;
Ad(4:6,1:3) = S*R;

end

%--------------------------------------------------------------------------
function y = isom(v)
y = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
end