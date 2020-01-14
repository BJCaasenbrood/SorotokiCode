function Background(Color)
if nargin < 1
Color = 'w';
end
fig = gcf;
set(fig,'color', Color);
fig.InvertHardcopy = 'off';
end

