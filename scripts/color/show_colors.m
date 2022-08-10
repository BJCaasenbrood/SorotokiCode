clr; 
%% loop over colormaps
figure(102);

for ii = 1:10
   subplot(10,1,ii) 
   showColormap(col(ii));
end

HTMLcode = col