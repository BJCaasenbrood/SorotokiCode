function plotpoint(p,varargin)
if length(p) == 3
plot3(p(1), p(2), p(3),varargin{:});
else
hold on; plot3(p(:,1), p(:,2),varargin{:});
end
end
