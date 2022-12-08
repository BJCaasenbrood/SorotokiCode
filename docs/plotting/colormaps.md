# Colormaps

##### Overview
``` matlab
fig(102,[10.75,10.75]);

colormaplist = {'blackwhite', 'bluesea', 'heatmap', 'inferno', 'metro', ...
                'noir','turbo','viridis','bounce','barney','evolution', ...
                'rainbow','polarmap','redgreen','soapbubble'};

N  = numel(colormaplist);
Nx = ceil(sqrt(N));
Ny = ceil(sqrt(N));

for ii = 1:N
   ax = subplot(Nx,Ny,ii) ;
   colorwheel();
   map = str2func(['@(x) ',colormaplist{ii}]);
   colormap(ax,map([]));
   axis off tight;
end
```

???+ "List of colormaps"
    <figure markdown>
    ![Image title](./img/colormaps.png){ width=800 }
    <figcaption>Improved colormaps of Sorotoki</figcaption>
    </figure>

##### Colormap manipulation