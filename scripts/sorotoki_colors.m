clr; 
%% loop over colors
figure(101);

x = linspace(0,2*pi,1e3);
y = @(x,a) sin(x + (a-1)/(pi));

for ii = 1:12
   plot(x,y(x,ii),'LineW',3); hold on;
   axis off;
   axis tight;
end


figure(102);

for ii = 1:12
   subplot(12,3,3*ii-2:3*ii-1) 
   showColormap(col(ii));
   subplot(12,3,3*ii) 
   text(0,0.33,['col(',num2str(ii),')'],'Fontsize',24);
   axis off;
   axis tight;
end
