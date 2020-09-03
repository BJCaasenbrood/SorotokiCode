function background(Color)
%BACKGROUND.m set the background color of the current figure;

if nargin < 1
    Color = 'w';
end
fig = gcf;
set(fig,'color', Color);
fig.InvertHardcopy = 'off';

end

