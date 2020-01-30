clc; clear; 

global Phi NModal NDof Ba Bc q

NModal = 1;
I6 = eye(6);
set = 1:6;
xa =  set(logical([0,1,1,1,0,0]));
xc = setdiff(set,xa);
Ba = I6(:,xa);
Bc = I6(:,xc);
Phi = matlabFunction(ShapeFunction(sym('t')));
NDof = size(Ba,2);

X = 0:0.01:1;
% q [t1,t2,t3,k4,k5,k6,e7,e8,e9]
q = zeros(NModal*NDof,1);
%q(3) = 3;
q(1) = 2; 
q(2) = 2; 
q(3) = .5;


% forward integration
g0 = [1,0,0,0,0,0,0];

[~,yf] = ode45(@(t,x) ForwardODE(t,x),X,g0);

hold on;
plot3(yf(:,7),yf(:,6),yf(:,5),'linewidth',1);
axis equal
view(30,30);
axis tight;

function dg = ForwardODE(t,g)
global Phi Ba q

ee = Ba*Phi(t)*q(:) + [0,0,0,1,0,0].';
Gamma = ee(4:6);
Kappa = ee(1:3);

Q = g(1:4);
R = Quat2Rot(Q);
A = StrainMap(R*Kappa(:));

dg = 0*g;
dg(1:4) = ((2*norm(Q))^(-1))*A*Q;
dg(5:7) = R*Gamma(:);
end

%---------------------------------------------------------------------- set
function A = StrainMap(K)
k1 = K(1); k2 = K(2); k3 = K(3);
A = [ 0, -k1, -k2, -k3; k1,   0, -k3,  k2; 
     k2,  k3,   0, -k1; k3, -k2,  k1,  0];
end

%---------------------------------------------------------------------- set
function R = Quat2Rot(q)
w = q(1); x = q(2); y = q(3); z = q(4);
Rxx = 1 - 2*(y^2 + z^2); Rxy = 2*(x*y - z*w); Rxz = 2*(x*z + y*w); 
Ryx = 2*(x*y + z*w); Ryy = 1 - 2*(x^2 + z^2); Ryz = 2*(y*z - x*w );
Rzx = 2*(x*z - y*w ); Rzy = 2*(y*z + x*w ); Rzz = 1 - 2 *(x^2 + y^2);

R = [Rxx, Rxy, Rxz; Ryx, Ryy, Ryz; Rzx, Rzy, Rzz];
end

%---------------------------------------------------------------------- set
function P = ShapeFunction(X)
global NModal NDof

for ii = 1:NModal, P(1,ii) = X.^(ii-1); end
for ii = 1:NDof, Pc{ii,1} = P; end

P = blkdiag(Pc{:}) + 1e-16*X;
end
