function c = cSquareSpiral(xc,yc,zc,r,d,N)   

R_ = (r/WaveFormSin(0.5*pi))*0.5;
Ri = R_ - d;
D = d/(2*pi*abs(N)); 
t = (2*pi*abs(N));

R = [cos(pi/4),-sin(pi/4), 0; sin(pi/4), cos(pi/4), 0; 0,0,1];
if N>0
    c = @(P) [(r - D*t*P).*WaveFormCos(t*P) + xc, (r - D*t*P).*WaveFormSin(t*P) + yc, P*0 + zc]*R;
else
    c = @(P) [(Ri + D*t*P).*WaveFormSin(t*P+0.5*pi) + xc, (Ri + D*t*P).*WaveFormCos(t*P+0.5*pi) + yc, P*0 + zc]*R;
end


end

function r = Radius(P,a,b,n)
X = abs(cos(P)/a).^n;
Y = abs(sin(P)/b).^n;
r = (X+Y).^(-1/n);
end

function y = WaveFormCos(x)
delta = 1e-12;
delta = delta*(1-x(:)) + 1e-12;
y = 1 - 2*acos((1-delta).*cos(x(:)))/pi;
end

function y = WaveFormSin(x)
delta = 1e-12;
delta = delta*(1-x(:)) + 1e-12;
y = 1 - 2*acos((1-delta).*sin(x(:)))/pi;
end



    