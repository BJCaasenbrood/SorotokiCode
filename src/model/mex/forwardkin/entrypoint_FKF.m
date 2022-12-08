L = 100; % length of robot
N = 120;  % number of discrete points on curve
M = 3;   % number of modes
X = linspace(0,L,N)';
Y = zeros(N,M);

for ii = 1:M
   Y(:,ii) = chebyshev(X/L,ii-1); 
end

Modes = [0,M,0,M,0,0];   % curvature-elongation

shp = Shapes(Y,Modes,'L0',L);

% shp.E    = 5.00;    % Young's modulus in Mpa
% shp.Nu   = 0.49;     % Poisson ratio
% shp.Rho  = 1000e-12; % Density in kg/mm^3
% shp.Zeta = 1;     % Damping coefficient
% shp = shp.rebuild();

%% eval theta
Fnc = @(x) shp.Theta(x);
s   = sort([shp.Sigma*shp.L0,...
            shp.Sigma*shp.L0+(2/3)*shp.ds]);

[nx,ny] = size(Fnc(0));
Theta   = zeros(nx,ny,numel(s));
Xi0     = zeros(6,1,numel(s));

for ii = 1:numel(s)
    Theta(:,:,ii) = Fnc(s(ii));
    Xi0(:,1,ii) = [0;0;0;1;0;0];
end


%% entries
x    = 0.01*rand(shp.NDim,1);
dx   = 0.01*rand(shp.NDim,1);
Ba   = shp.Ba;
ds   = shp.ds;
p0   = zeros(3,1);
Phi0 = eye(3);
xia0 = Xi0;
Th   = Theta;

[g,J] = computeForwardKinematicsFast(x,dx,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zeroclc
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ba);

%clearvars -except g J