clr; 
%% loop over colormaps
figure(102);

wheel();
colormap(viridis(0));

% for ii = 1:12
%    subplot(12,1,ii) 
%    showColormap(col(ii));
% end
% 
% HTMLcode = col


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