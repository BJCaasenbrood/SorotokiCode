clr;
%

N = 12; X = linspace(0,pi*3,1e5);
hold off;
for ii=1:N
    Y = cos(X+0.75*ii*pi/N);
    plot(X,Y,'color',col(ii),'linewidth',3);
    hold on;
end