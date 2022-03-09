function y = pdac(x,varargin)
    
    if isempty(varargin)
       beta = 0;
    else
       beta = varargin{1};
    end
    
    y = clamp((x+100+beta)/200,0,1);
end

