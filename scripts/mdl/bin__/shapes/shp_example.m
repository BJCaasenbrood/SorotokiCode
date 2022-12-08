clr; beep off;
%% settings
L = 1;   % manipulator length
M = 6;   % number of modes
N = 500; % grid on SR

%% build basis
%x   = linspace(0,1,N).';                % 
%Y   = GenerateFunctionSpace(x,N,M,L);   % 

Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,0,0,0,0]);          % generate basis

q = zeros(M,1);
q(1) = -1;
q(2) = -3;
q(4) = 3;

%% render curve
p = shp.FK(q);

sdf = sPolyline(p(:,[1,3]));
sdf.BdBox = [0,1,-0.5,0.2];
sdf.show('Quality',300);
colormap(viridis(0))

hold on;

fplot(p(end,[1,3]),'w.','MarkerSize',25);
axis equal;
