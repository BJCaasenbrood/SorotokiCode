function y = smoothfall(x)

I1 = (x<=0);
I2 = (x>0 & x<=1);
I3 = ~I1 & ~I2;

y = x;
y(I1) = 1;
y(I2) = -3*x(I2).^2 +2*x(I2).^3 + 1;
y(I3) = 0;

end