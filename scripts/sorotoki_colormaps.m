%% SOROTOKI: colormap overview
clr; 
%% loop over colormaps
fig(102,[10.75,10.75]);

name = {'blackwhite', 'bluesea', 'heatmap', 'inferno', 'metro', ...
        'noir','turbo','viridis','bounce','barney','evolution','rainbow',...
        'polarmap','redgreen','soapbubble'};

N  = numel(name);
Nx = ceil(sqrt(N));
Ny = ceil(sqrt(N));

for ii = 1:N
   ax = subplot(Nx,Ny,ii) ;
   colorwheel();
   map = str2func(['@(x) ',name{ii}]);
   colormap(ax,map(1));
   axis off;
   axis tight;
   title(['\texttt{',name{ii},'}'],'Fontsize',14);
end

