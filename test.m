clr;

[x,y,z] = peaks(300);

z = z*0.25;

figure
surf(x,y,z,'linestyle','none');

colormap(viridis)

axis tight
axis equal
axis off;
