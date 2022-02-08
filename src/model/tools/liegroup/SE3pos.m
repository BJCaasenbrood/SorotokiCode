function p = SE3pos(X,varargin)

if size(X,1) == 7
    q = R7(1:4,1);     % quaterions
    p = R7(5:7,1);     % position
elseif numel(X) == 9
    if ~isempty(varargin)
       p = varargin{1};
    else
       p = zeros(3,1); 
    end
end

end

