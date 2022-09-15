clr; beep off;
%% settings
L = 1;   % manipulator length
M = 3;   % number of modes
N = 500; % grid on SR

%% build basis
x   = linspace(0,1,N).';                % 
Y   = GenerateFunctionSpace(x,N,M,L);   % 
shp = Shapes(Y,[M,M,0,0,0,0]);          % generate basis

q = zeros(2*M,1);
q(1) = 17;
q(5) = 15;

%% render curve
p = shp.FK(q);
[~,Et] = shp.tangentPoint(q);

subplot(2,1,1);
cplot3(p(:,1),p(:,2),p(:,3),Et,turbo,'LineW',12);
axis equal;

subplot(2,1,2);
plot(x,Et,'LineW',1.5);

function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
   %Y(:,ii) = pcc(X/L,ii,M); %chebyshev(X/L,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end
