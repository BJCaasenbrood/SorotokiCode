function c = cHelix(xc,yc,zc,r,l,N)    
  
c = (l/(2*pi*N));
t = ((2*pi)*N);
    
c = @(P) [r*cos(t*P) + xc, r*sin(t*P) + yc, c*t*P + zc];


    