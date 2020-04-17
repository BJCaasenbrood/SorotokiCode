function y = legendre(x,n)
x = (2*x-1);

if (n==0),      y = 1 + x*0;
elseif(n==1),   y = x;
elseif(n==2),   y = 0.5*(3*x.^2 - 1);
elseif(n==3),   y = 0.5*(5*x.^3 - 3*x);
else,           y = 1 + x*0;
end
end


