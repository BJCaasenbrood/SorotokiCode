N = 10; 
X = linspace(0,pi*3,1e5);
label = [];

for ii=1:N
    Y = cos(X+1*ii*pi/N);
    plot(X,Y);
    hold on;
    pause(0.01);
    label{ii} = ['col(',num2str(ii),')'];
end

% magick color change
sorocolor;

xlim([0,14]);
legend(label);
axis off;