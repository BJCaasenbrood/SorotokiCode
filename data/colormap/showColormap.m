function showColormap(x)
if nargin < 1, x = viridis; end

[X,Y] = meshgrid(1:500,1:40);

surf(X,Y,X,'linestyle','none'); hold on;
plot3([1,500,500,1,1],[1,1,40,40,1],...
       600*ones(1,5),'k','linewidth',0.1);
   
colormap(gca,x); view(0,90);
axis tight; axis off;
background(gitpage);

end

