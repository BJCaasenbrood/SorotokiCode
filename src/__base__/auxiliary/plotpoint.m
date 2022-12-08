function plotpoint(p,varargin)
% PLOTPOINT.m plots point or pointcloud in 3D space. Usage:
% 
MarkerStyle = '.';
if numel(p) == 3
    plot3(p(1), p(2), p(3),'Marker',MarkerStyle, ...
        'MarkerSize',35,varargin{:});
else
    hold on; plot3(p(:, 1), p(:, 2), varargin{:});
end
end
