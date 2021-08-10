function  Colors = colorbase
I = imread('colorwheel.PNG');
W = size(I,1);
R = 0.45*W;
Dth = 2*pi/12;
% image(I); hold on;
% for i = 1:12
%    plot(W/2+R*cos(xth(i)),W/2 + R*sin(xth(i)),'w.','markersize',10);
% end
Colors = {};
loopset = cumsum([0,4,6,8]);
Id0 = 1; % light blue color
empset = [];
for i = 1:2
    th = Dth*(Id0+loopset);
    empset = [empset,Id0+loopset];
    X = W/2+R*cos(th);
    Y = W/2+R*sin(th);
    %plot(X,Y,'w.-','markersize',10);
    Id0 = Id0+3;
    for j = 1:4
       py = floor(X(j));
       px = floor(Y(j));
       Colors{end+1} = rgb2hex([I(px,py,1),I(px,py,2),I(px,py,3)]);
    end
end

nonset = setdiff(1:12,mod(empset,12));

for j = nonset
    th = Dth*j;
    X = W/2+R*cos(th);
    Y = W/2+R*sin(th);
    py = floor(X);
    px = floor(Y);
    Colors{end+1} = rgb2hex([I(px,py,1),I(px,py,2),I(px,py,3)]);
end

end

