
function c = cCircle(P,xc,yc,zc,r)   
  
t = (2*pi)*P;

c = @(P) [r*cos(t) + xc, r*sin(t) + yc, zc];


    