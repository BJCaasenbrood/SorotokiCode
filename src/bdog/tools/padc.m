function y = padc(x,varargin)
    
    if isempty(varargin)
       beta = 0;
    else
       beta = varargin{1};
    end
    
    y = (x - 5)*20 - beta;
end

