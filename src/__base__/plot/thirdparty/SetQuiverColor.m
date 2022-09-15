function SetQuiverColor(q,currentColormap,varargin)
%--------------------------------------------------
% function SetQuiverColor(q,currentColormap)
%
% INPUT:
%   q = handle to quiver plot
%   currentColormap = e.g. jet;
% OPTIONAL INPUT ('Field',value):
%   'range' = [min,max]; % Range of the magnitude in the colorbar
%                          (used to possibly saturate or expand the color used compared to the vectors)
%   'mags' = magnitude; % Actual magnitude of the vectors
%
% Example:
%   [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%   z = x .* exp(-x.^2 - y.^2);
%   [u,v,w] = surfnorm(x,y,z);
%   q = quiver3(x,y,z,u,v,w);
%   mag = 1+3.*rand(size(u));   % Creates number between 1 and 4
%   colormap(jet);
%   colorbar;
%   SetQuiverColor(q,jet,'mags',mag,'range',[-2 8]);  % Color range between -2 8 => all colors are not used
%   caxis([-2 8]);
%   set(gca,'Color','k');
%
%--------------------------------------------------
%   Authorship:
%     This code is heavily based from the answer by the user Suever on Stackoverflow forum
%     at: https://stackoverflow.com/questions/29632430/quiver3-arrow-color-corresponding-to-magnitude
%
%     I, Alexandre De Spiegeleer, only added minor changes to the original answer to have a bit more flexibility.
%--------------------------------------------------

%// Set default values
range = [];
mags = [];

%// Read the optional range value
if find(strcmp('range',varargin))
  range = varargin{ find(strcmp('range',varargin))+1 };
end

qU = q.UData(~isnan(q.UData));
qV = q.VData(~isnan(q.VData));
qW = q.WData(~isnan(q.WData));

%// Compute/read the magnitude of the vectors
if find(strcmp('mags',varargin))
  mags = varargin{ find(strcmp('mags',varargin))+1 };
  mags = mags(~isnan(mags)&~isnan(q.UData));  % This reshapes automatically
else
  mags = sqrt(sum(cat(2, qU, qV, ...
             reshape(qW, numel(qU), [])).^2, 2));
end
%// If range is auto, take range as the min and max of mags
if isstr(range) & strcmp(range,'auto')
  range = [min(mags) max(mags)];
end

%// Change value depending on the desired range
if ~isempty(range) & isnumeric(range) & numel(range)==2
  range = sort(range);
  mags(mags>range(2)) = range(2);
  mags(mags<range(1)) = range(1);
end

%// Now determine the color to make each arrow using a colormap
if ~isempty(range) & isnumeric(range) & numel(range)==2
  Edges = linspace(range(1),range(2),size(currentColormap, 1)+1);
  [~, ~, ind] = histcounts(mags, Edges);
else
  [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
end

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// Color data
cd_head = reshape(cmap(1:3,:,:), [], 4).';
cd_tail = reshape(cmap(1:2,:,:), [], 4).';

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', cd_head);

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', cd_tail);
