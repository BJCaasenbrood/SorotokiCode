function y = smoothfall(x)

I1 = (x<=0);
I2 = (x>0 & x<=1);
I3 = ~I1 & ~I2;

y = x;
y(I1) = 1;
y(I2) = -3*x(I2).^2 +2*x(I2).^3 + 1;
y(I3) = 0;

% for ii = 1:length(X)
%     x = X(ii);
%     if x<=0, y(ii,1) = 0;
%     elseif (x>0 && x<=1), y(ii,1) = 3*x.^2 -2*x.^3;
%     else, y(ii,1) = 1;
%     end
% end
end