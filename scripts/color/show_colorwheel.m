clr; 
%% loop over colormaps
figure(102);

wheel();
colormap(cmap_orange);

%% color wheel

function [X,Y,C] = wheel
n = 1e3;
m = round(0.2*n);
r = (m:n)'/n;
theta = pi*(-n:n)/n;
X = r*cos(theta);
Y = r*sin(theta);
C = atan2(Y,X);


s = pcolor(X,Y,C);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
axis equal;
axis off;
colorbar;
end