%% SOROTOKI: colormap overview
clr; 
%% loop over colormaps
figure(102);

name = {'blackwhite', 'bluesea', 'heatmap', 'inferno', 'metro', ...
    'noir','turbo','viridis','bounce','barney','evolution','rainbow',...
    'polarmap','redgreen','soapbubble'};

N  = numel(name);
Nx = ceil(sqrt(N));
Ny = ceil(sqrt(N));

for ii = 1:N
   ax = subplot(Nx,Ny,ii) ;
   %showColormap(map{ii});
   colorwheel();
   map = str2func(['@(x) ',name{ii}]);
   colormap(ax,map(1));
   axis off;
   axis tight;
   title(name{ii});
end

background(gitpage)
