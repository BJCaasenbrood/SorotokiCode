clr; 
%% loop over colormaps
figure(101);

cm = @(x) turbo(x);
map = {cm([]),cm(5),cm(25),cm(500),cm(0),cm(-1)};

opts = {'','5','25','500','0','-1'};

for ii = 1:length(opts)
   subplot(length(opts),3,3*ii-2:3*ii-1) 
   showColormap(map{ii});
   subplot(length(opts),3,3*ii) 
   text(0,0.33,['viridis(',opts{ii},')'],'interpreter','tex');
   axis off;
   axis tight;
end

background(gitpage)
