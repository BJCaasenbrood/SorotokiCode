clr; 
%% loop over colormaps
figure(102);

for ii = 1:12
   subplot(12,3,3*ii-2:3*ii-1) 
   showColormap(col(ii));
   subplot(12,3,3*ii) 
   text(0,0.33,['col(',num2str(ii),')'],'interpreter','tex');
   axis off;
   axis tight;
end
