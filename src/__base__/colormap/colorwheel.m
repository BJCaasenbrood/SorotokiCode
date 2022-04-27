function colorwheel()

n = 60;
r = (10:n)'/n;
theta = pi*(-n:n)/n - pi/2;
X = r*cos(theta);
Y = r*sin(theta);
C = (0*r + 1)*theta;

h = pcolor(X,Y,C);
h.LineStyle = 'none';
axis equal tight;
axis off;

end

