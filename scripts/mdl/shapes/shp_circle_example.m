clr; beep off;
%% settings
L = 1;   % manipulator length
M = 3;   % number of modes
N = 301; % grid on SR

%% build basis
x   = linspace(0,1,N).';                % 
Y   = GenerateFunctionSpace(x,N,M,L);   % 
shp = Shapes(Y,[0,M,0,0,0,0]);          % generate basis

shp.Xi0 = @(x) [0;2*pi;0;1;0;0];
shp = shp.rebuild();

%%
q = zeros(M,1);
q(1) = 0;

%% render curve
p = shp.FK(q);
[~,Et] = shp.tangentPoint(q);

subplot(2,1,1);
cplot3(p(:,1),p(:,2),p(:,3),Et,turbo,'LineW',5);
axis equal;

subplot(2,1,2);
plot(x,Et,'LineW',1.5);


function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);
jj = 1;
for ii = 2:2:2*M
   Y(:,jj) = chebyshev(X/L,jj); % chebyshev
   %Y(:,ii) = pcc(X/L,ii,M); %chebyshev(X/L,ii-1); % chebyshev
   jj = jj + 1;
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end
