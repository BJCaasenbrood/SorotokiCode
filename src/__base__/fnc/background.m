function background(Color)
% function background(color)
%--------------------------------------------------------------------------
% Function set the background color of the current figure (i.e. gcf);
%
% Usage:
% background([1,1,1]); % set background to white
% background('w');     % set background to white
%--------------------------------------------------------------------------

if nargin < 1
    Color = 'w';
end

fig = gcf;
set(fig,'color', Color);
fig.InvertHardcopy = 'off';

end

