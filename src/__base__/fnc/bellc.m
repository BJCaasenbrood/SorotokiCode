function y = bellc(x,varargin)
if ~isempty(varargin)
    if numel(varargin) == 1
        sigma = varargin{1};
        mu = 0;
    elseif numel(varargin) == 2
        sigma = varargin{1};
        mu    = varargin{2};
    end
else
    sigma = 1;
    mu = 0;
end
    
beta = -1/(2*sigma^2);
y = (1/(sigma*2*pi))*exp(beta*(x-mu).^2);

end

