% %a = [2.23e3,1.74e3,-4.55e2];
% a = [4.23e-1,3.99e-1,-2.29e-1];
% 
% x = linspace(-0.8,1.2,1e3);
% 
% plot(x,hyp(x,a),'linewidth',1.5)
% 
% 
% function y = hyp(x,a)
% y = a(1) + a(2)*(tanh(a(3)*x).^2 - 1);
% end

xq = linspace(x(1),x(end),50);
vq = interp1(x,y,xq);

plot(0.35+xq(2:end),diff(vq)./diff(xq),'b.')
