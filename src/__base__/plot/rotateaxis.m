function rotateaxis(varargin)
if isempty(varargin)
   deg = 90;
else
   deg = varargin{1};
end
AxesH    = gca;
UpVector = [-sind(deg), cosd(deg), 0];
DAR      = get(AxesH, 'DataAspectRatio');
set(AxesH, 'CameraUpVector', DAR .* UpVector);
end

