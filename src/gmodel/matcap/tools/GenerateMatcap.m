function GenerateMatcap(name,cmap)
close all;
path = cdsoro;
pathmcap = [path,'\data\matcap\'];
pathimg = [path,'\data\matcap\'];

[x,y,z] = sphere(100);
z = (z);

Dx = 512;
Dy = 512;
b = 0.5715;
a = 512;
figure('Position', [a,a,b*Dx,b*Dy]);
cplane(x,y,z.^2);
%image(unique(x(:)),flipud(unique(y(:))),z*255);
axis equal;
axis([-1 1 -1 1]);
ax = gca;
%AxesH = axes;
ax.Clipping = 'off';
%InSet = get(AxesH, 'TightInset');
%set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
%set(gcf, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
axis image
view(0,90);
caxis([0 1]);
shading interp;
colormap(cmap);
background(0.85*cmap(1,:));
axis off;
drawnow;
F = getframe;
imwrite(F.cdata,name);
end

