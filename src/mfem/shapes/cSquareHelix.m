function c = cSquareHelix(xc,yc,zc,a,b,n,l,N)    
c = (l/(2*pi*N));
t = ((2*pi)*N);
c = @(P) [Radius(t*P,a,b,n).*cos(t*P) + xc, Radius(t*P,a,b,n).*sin(t*P) + yc, c*t*P + zc];
%c = @(P) [r*cos(t*P) + xc, r*sin(t*P) + yc, c*t*P + zc];
end

function r = Radius(P,a,b,n)
X = abs(cos(P)/a).^n;
Y = abs(sin(P)/b).^n;
r = (X+Y).^(-1/n);
end