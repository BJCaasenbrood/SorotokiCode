clr; 
%% loop over colormaps
figure(101);

map = {blackwhite, bluesea, heatmap, inferno, metro, ...
       noir,turbo,viridis,bounce,barney};

name = {'blackwhite', 'bluesea', 'heatmap', 'inferno', 'metro', ...
        'noir','turbo','viridis','bounce','barney'};

for ii = 1:length(map)
   subplot(length(map),3,3*ii-2:3*ii-1) 
   showColormap(map{ii});
   subplot(length(map),3,3*ii) 
   text(0,0.33,name{ii},'interpreter','tex');
   axis off;
   axis tight;
end

background(gitpage*0.75)
