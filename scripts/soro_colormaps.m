clr; 
%% loop over colormaps
figure(102);

map = {blackwhite, bluesea, heatmap, inferno, metro, ...
    noir,turbo,viridis};

for ii = 1:length(map)
   subplot(length(map),1,ii) 
   showColormap(map{ii});
   axis tight;
end
