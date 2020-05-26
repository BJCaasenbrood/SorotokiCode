clr; 
%% loop over colormaps
figure(102);

for ii = 1:12
   subplot(12,1,ii) 
   showColormap(col(ii));
end
