clr;
x = linspace(-pi,pi,1e3);

%
fig(101,[9,9]);
plot(x,sin(x.^2));

xlabel('$x$');
ylabel('$sin(x^2)$');
xlim([-pi,pi]);
ylim([-2,2])