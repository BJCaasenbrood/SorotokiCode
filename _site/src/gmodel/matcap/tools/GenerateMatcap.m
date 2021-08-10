function GenerateMatcap(name,cmap)
path = cdsoro;
pathmcap = [path,'\data\matcap\'];
pathimg = [path,'\data\matcap\'];

[x,y,z] = sphere(100);
z = (z);

figure('Position', [100 100 500 500]);
surf(x,y,z.^2,'linestyle','none');
%image(unique(x(:)),flipud(unique(y(:))),z*255);
axis square;
axis([-1 1 -1 1]);
ax = gca;
ax.Clipping = 'off';
AxesH = axes;
drawnow;
% InSet = get(AxesH, 'TightInset');
%set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
 set(gcf, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
%axis image
view(0,90);
caxis([0 1]);
axis off;
shading interp;
colormap(cmap);
background(cmap(1,:));
axis off;
print(gcf,name,'-dpng','-r300');
end

