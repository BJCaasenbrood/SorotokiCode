[shp] = SetupModel();

p = shp.FK([2e-2,-0.1]);

figure(102); %clf;
plot(p(:,1),p(:,3),'bo-'); hold on;
axis equal; 

function [shp] = SetupModel()
L = 86;  % intrinsic length in mm
M = 1;   % number of modes
N = 15;  % grid on SR
x = linspace(0,L,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = pcc(x/L,ii,M);       % piece-wise constant
end

% soft sobotics shapes
shp = Shapes(Y,[0,M,0,M,0,0],'L0',L,'g0',SE3(roty(pi/2),zeros(3,1)));
%shp = shp.rebuild();

end