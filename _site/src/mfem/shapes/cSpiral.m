
function c = cSpiral(xc,yc,zc,r,d,N)   
Ri = r - d;
D = d/(2*pi*abs(N)); 
t = (2*pi*abs(N));

if N>0
    c = @(P) [(r - D*t*P).*cos(t*P) + xc, (r - D*t*P).*sin(t*P) + yc, P*0 + zc];
else
    c = @(P) [(Ri + D*t*P).*sin(t*P+0.5*pi) + xc, (Ri + D*t*P).*cos(t*P+0.5*pi) + yc, P*0 + zc];
end

    