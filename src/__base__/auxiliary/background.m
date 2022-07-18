function background(x)
% BACKGROUND overwrites the background color of the current figure, i.e.,
%   gcf. Altenatively, it can be used to assign a background image. Usage:
%
%   % setting colors
%   background([1,1,1]);    % set background to white
%   background('w');        % set background to white
%
%   % settting images as background
%   I = imread('egg.jpg');
%   background(I);          % makes egg.jpg the background
%
%   Last edit: July 15, 2022.

if nargin < 1
    x = 'w';
end

if (isa(x,'double') && size(x(:),1) == 3) || (isa(x,'char'))
    fig = gcf;
    set(fig,'color', x);
    fig.InvertHardcopy = 'off';
else
    ah = axes('unit', 'normalized', 'position', [0 0 1 1]);
    imagesc(x);
    set(ah,'handlevisibility','off','visible','off');
    uistack(ah,'bottom');
end

end


