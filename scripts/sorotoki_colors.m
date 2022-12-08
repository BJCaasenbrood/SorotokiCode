%clr; 
%% loop over colors
figure(102);

%colorlist = redgreen(10);
colorlist = col;
%%
x = linspace(0,2*pi,1e3);
y = @(x,a) sin(x + (a-1)/(.444*pi));
Leg = {};

for ii = 1:size(colorlist,1)
   plot(x,y(x,ii),'LineW',5,...
       'Color',colorlist(ii,:)); hold on;
   axis off;
   axis tight;
   Leg{ii} = num2str(ii);
end

legend(Leg,'Orientation','Vertical','Location','west');