clr;

fig(101,[11,9]);
subplot(2,5,[1:5]);
showcolor(color_primary,color_primary_darker,color_primary_darkest,color_base);
legend({'primary','primary darker','primary darkest','base'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

subplot(2,5,[6:10]);
showcolor(color_gray_warm_dark,color_gray_warm_light,color_primary_darkest,color_gray_cool_light);
legend({'gray warm dark','gray warm light','primary darkest','gray cool light'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{1} = imresize(f,[NaN,1024]);

clf;
fig(101,[10,7.5]);
%subplot(2,5,[1:5]);
showcolor(color_primary_alt_darkest,color_primary_alt_dark,...
    color_primary_alt,color_primary_alt_light,color_primary_alt_lightest);
legend({'... darkest','... dark','primary alt','... light',...
    '... lightest'},'Location','southoutside','Edgecolor','none',...
    'orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{2} = imresize(f,[NaN,1024]);

clf;
fig(101,[10,7.5]);
%subplot(2,5,[1:5]);
showcolor(color_secondary_darkest,color_secondary_dark,...
    color_secondary,color_secondary_light,color_secondary_lightest);
legend({'... darkest','... dark','secondary','... light',...
    '... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{3} = imresize(f,[NaN,1024]);

clf;
fig(101,[10,7.5]);
%subplot(2,5,[1:5]);
showcolor(color_gray_dark,color_gray,color_gray_light,...
    color_gray_lighter,color_gray_lightest);
legend({'... dark','gray','... light','... lighter','... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{4} = imresize(f,[NaN,1024]);

close all;
figure(101);
%out = imtile(F,'GridSize',[3,1]);
out = cat(1,F{:});
imshow(out);

F = [];
fig(102,[10,7.5]);
%subplot(2,5,[1:5]);
showcolor(color_gold,color_gold_light,color_gold_lighter,color_gold_lightest);
legend({'gold','... light','... lighter','... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{1} = imresize(f,[NaN,1024]);

clf
fig(102,[10,7.5]);
%subplot(2,5,[1:5]);
showcolor(color_green,color_green_light,color_green_lighter,color_green_lightest);
legend({'green','... light','... lighter','... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{2} = imresize(f,[NaN,1024]);

clf
fig(102,[10,7.5]);
%subplot(2,5,[1:5]);
showcolor(color_cool_blue,color_cool_blue_light,color_cool_blue_lighter,...
    color_cool_blue_lightest); 
legend({'... blue','... light','... lighter','... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{3} = imresize(f,[NaN,1024]);

clf
fig(102,[10,7.5]);
%subplot(2,5,[1:5]);
showcolor(color_visited,color_visited_light,color_visited_lighter,color_visited_lightest);
legend({'visited','... light','... lighter','... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{4} = imresize(f,[NaN,1024]);

clf
fig(102,[10,7.5]);
showcolor(color_orange,color_orange_light,color_orange_lighter,...
    color_orange_lightest);
legend({'orange','... light','... lighter','... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{5} = imresize(f,[NaN,1024]);

clf
fig(102,[10,7.5]);
showcolor(color_visited_alt,color_visited_alt_light,color_visited_alt_lighter,...
    color_visited_alt_lightest);
legend({'visited alt','... light','... lighter','... lightest'},...
    'Location','southoutside','Edgecolor','none','orientation','Horizontal');

pause(.2);
f = getframe(gcf).cdata;
F{6} = imresize(f,[NaN,1024]);

close(102);
figure(102);
%out = imtile(F,'GridSize',[3,1]);
out = cat(1,F{:});
imshow(out);