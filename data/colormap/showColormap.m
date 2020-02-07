function showColormap(x)
if nargin < 1, x = viridis; end

f = figure(221); clf;
%set(f, 'Position',  [100, 100, 500, 140]);

[X,Y] = meshgrid(1:500,1:40);

surf(X,Y,X,'linestyle','none'); hold on;
plot3([1,500,500,1,1],[1,1,40,40,1],600*ones(1,5),'w','linewidth',1)
colormap(x);
view(0,90);
axis equal; axis tight; axis off;
text(-5,-20,'0','Fontsize',14,'Color','w');
text(495,-20,'1','Fontsize',14,'Color','w');
background('k')

end

