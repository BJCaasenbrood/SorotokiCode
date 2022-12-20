%#codegen
function [M,C,K,R,G,B,p,Phi,J,dJdt,Vg,Kin] = computeLagrangianFast(x,dx,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zero
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ms,...      % evaluated muscle vector
    Ba,...      % state to strain matrix
    Ktt,...     % geometric stiffness
    Mtt,...     % geometric inertia
    Zeta,...    % dampings coefficient
    Gvec)       % gravitional vector  

% compute total length
n = numel(x);
m = size(Ms,3);
s = 0;

Z1 = zeros(6,6+2*(n-1));
Z2 = zeros(n,3*n+1);
Z3 = zeros(n,m);

Z1(1:3,1:3) = Phi0;
Z1(1:3,4)   = p0;

%NLStiff = false; 

for ii = 1:(size(Th,3)/2)
    
    % first EL-diff eval
    [K1Z1,K1Z2,K1Z3] = LagrangianODEX(x, dx, Z1,...
        Th(:,:,2*ii-1), xia0(:,1,2*ii-1), Ms(2*ii-1,:,:), ...
        Ba, Mtt, Ktt, Gvec);
    
    % second EL-diff eval
    [K2Z1,K2Z2,K2Z3] = LagrangianODEX(x, dx, Z1 + (2/3)*ds*K1Z1,...
        Th(:,:,2*ii), xia0(:,1,2*ii), Ms(2*ii,:,:), ...
        Ba, Mtt, Ktt, Gvec);
    
    % update integrands
    s  = s  + ds;
    Z1 = Z1 + 0.25*ds*(K1Z1 + 3*K2Z1);
    Z2 = Z2 + 0.25*ds*(K1Z2 + 3*K2Z2);
    Z3 = Z3 + 0.25*ds*(K1Z3 + 3*K2Z3);

end

% recover the kinematics entities
p    = Z1(1:3,4);
Phi  = Z1(1:3,1:3);
B1   = Z1(1:6,5:5+n-1);
B2   = Z1(1:6,6+n-1:6+2*(n-1));
J0   = Admapinv(Phi,p)*B1;
dJ0  = Admapinv(Phi,p)*B2;

V    = J0*dx; 
adV  = admap(V);

% set linear velocity zero -- dJ worldframe 
V(4:6)  = 0;
adV_ = admap(V);
J    = Admap(Phi,[0;0;0])*J0;           % world-frame Jacobian
dJdt = Admap(Phi,[0;0;0])*adV_*J0 + ... % world-frame dJacobian/dt
       Admap(Phi,[0;0;0])*(-adV*J0 + dJ0);

% recover the dynamics entities
M  = Z2(1:n,1:n);
C  = Z2(1:n,n+1:2*n);
K  = Z2(1:n,2*n+1:3*n);
G  = Z2(1:n,3*n+1);

Vg  = Z1(5,4);
Kin = Z1(6,4);

R = Zeta*K;
B = Z3;
end

function [dZ1,dZ2,dZ3] = LagrangianODEX(x,dx,Z1,...
    Theta,xia0,Mscl,Ba,Mtt,Ktt,Gvec)

n     = numel(x);
m     = size(Mscl,3);
p_    = Z1(1:3,4);
Phi_  = Z1(1:3,1:3);
J_    = Z1(1:6,5:5+n-1);
Jt_   = Z1(1:6,6+n-1:6+2*(n-1));

%build geometric strain
XI = Ba*Theta*x + xia0;

Gamma = XI(1:3);
U     = XI(4:6);

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*isomSO3(Gamma);

A   = Admap(Phi_,p_);
Ai  = Admapinv(Phi_,p_);

% build jacobian
Jg   = Ai*J_;
Jgt  = Ai*Jt_;
V    = Jg*dx;
adV  = admap(V);

BTh = Ba*Theta;

dJ  = A*BTh;
dJt = A*adV*BTh;

% compute inertia, coriolis, gravity
dM = (Jg).'*Mtt*Jg;
dC = (Jg).'*((Mtt*adV - adV.'*Mtt)*Jg  + Mtt*Jgt);
dG = -(Jg).'*(Ai*Mtt*[0;0;0;Gvec]);

% compute (nonlinear stiffness)
dK = (((BTh).'*Ktt*(BTh)));

% compute grav. potential energy
dKe = 0.5*V.'*Mtt*V;
dVg = Mtt(4,4)*p_.'*Gvec;

dZ1                      = zeros(6,6+2*(n-1));
dZ1(1:3,1:3)             = dPhi;
dZ1(1:3,4)               = dp;

dZ1(1:6,5:5+n-1)         = dJ;
dZ1(1:6,6+n-1:6+2*(n-1)) = dJt;

dZ1(5,4)                 = dVg;
dZ1(6,4)                 = dKe;

dZ2 = zeros(n,3*n+1);
dZ2(1:n,1:n)             = dM;
dZ2(1:n,n+1:2*n)         = dC;
dZ2(1:n,2*n+1:3*n)       = dK;
dZ2(1:n,3*n+1)           = dG;

dZ3 = zeros(n,m);

for ii = 1:m
    
    di  = [Mscl(:,1:3,ii).';1]; % lifted position vector
    ddi = [Mscl(:,4:6,ii).';1]; % lifted position vector deriv.
    
    XiHat      = hat(XI);
    tangentRaw = XiHat*di + ddi;
    ti         = tangentRaw(1:3)/norm(tangentRaw(1:3));
    
    Fu        = [hat(di(1:3))*ti;ti];
    dZ3(:,ii) = ((BTh.')*Fu);
end

end
%--------------------------------------------------------------------------
function y = isomSO3(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; 
     x3, 0, -x1; 
     -x2, x1, 0];
end
%--------------------------------------------------------------------------
function Y = hat(X)

if numel(X) == 6
    S = isomSO3(X(1:3));
    U = X(4:6);
    
    Y = zeros(4,4);
    Y(1:3,1:3) = S;
    Y(1:3,4)   = U;
    
elseif numel(X) == 3
    Y = isomSO3(X);
end
    
end

%--------------------------------------------------------------------------
function g = admap(x)
    W = x(1:3);
    U = x(4:6); 
    g = zeros(6);
    Wh = isomSO3(W); 
    Uh = isomSO3(U);
    g(1:3,1:3) = Wh;
    g(4:6,4:6) = Wh;
    g(4:6,1:3) = Uh;
end
%--------------------------------------------------------------------------
function Ad = Admap(R,p)
    S = isomSO3(p);
    Ad = zeros(6);
    Ad(1:3,1:3) = R;
    Ad(4:6,4:6) = R;
    Ad(4:6,1:3) = S*R;
end
%--------------------------------------------------------------------------
function Ad = Admapinv(R,p)
    Rt = R.';   
    S = isomSO3(p);
    Ad = zeros(6);
    Ad(1:3,1:3) = Rt;
    Ad(4:6,4:6) = Rt;
    Ad(4:6,1:3) = Rt*S.';
end