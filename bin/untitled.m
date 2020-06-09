clr
%%
x = linspace(0+1e-6,1-1e-6,1e3);
% 
for ii = 1:3
y1 = sign(z1(x)).*chebyshev(z1(x),ii-1);
y2 = sign(k2(x)).*chebyshev(k2(x),ii-1);

hold on;
subplot(3,2,2*ii-1);
shade(x,y1,'Color',col(1),'Linewidth',1.5);
axis([0 1 -1 1]); grid on;
set(gca,'linewidth',0.75)
subplot(3,2,2*ii);
shade(x,y2,'Color',col(2),'Linewidth',1.5);
axis([0 1 -1 1]); grid on;
set(gca,'linewidth',0.75)
end

% plot(x,k1(x))

function y = z1(x)
%y = x;
y(x>0.5) = 0;
y(x<=0.5) = 2*x(x<=0.5);
end

function y = z2(x)
%y = x-1;
y(x<0.5) = 0;
y(x>=0.5) = 2*x(x>=0.5);
end

function y = k1(x)
z = max(sign(-2*x+1),0);
y = 2*z.*x;
end

function y = k2(x)
y = max(2*(x)-1,0);
end

function plota(x,y,n)
plot(x,y,'color',col(n),'linewidth',1.5); hold on;
area(x,y,'facealpha',0.15,'facecolor',col(n),'linestyle','none');
end