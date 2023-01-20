%CPLOT - Plots data points with color map
%
%   h = cplot(x, y, c, map, varargin)  plots the data points with coordinates (x, y) 
%   with colors specified by the c vector, and color map defined by the input map.
%   The option varargin allows additional plotting parameters to be passed
%   to the standard plot routine within Matlab.
%
%   Output h is a list of plot handles
%
%   Example:
%       x = linspace(0, 6*pi,100);
%       y = sin(x);
%       c = cos(x);
%       map = turbo;
%       h = cplot(x, y, c, map, 'LineWidth', 2);
%
%   Note:
%      Vectors X, Y and C must be the same size
%
%   See also CPLOT3, CPLOTARROW

function h = cplot(x,y,c,map,varargin)

if ~(all(size(x) == size(y)) && all(size(x) == size(c)))
    error('Vectors X,Y and C must be the same size');
end

N = size(map,1);
cmax = max(c);
cmin = min(c);
cint = (cmax-cmin)/N;

indices = 1:length(c);
status  = ishold;               

for k = 1:N
    ii = logical(c >= cmin+k*cint) + logical(c <= cmin+(k-1)*cint);
    jj = ones(size(ii)); jj(1:end-1) = ii(2:end);
    ii = logical(ii .* jj);
    X = x; X(indices(ii)) = NaN;
    Y = y; Y(indices(ii)) = NaN;
    h(k) = plot(X,Y,'Color',map(k,:),varargin{:});
    hold on;
end

if status == 0
    hold off; 
end   

end

