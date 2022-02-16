%#codegen
function [M,C,K,R,G,p,Phi,J,Vg,Kin] = computeLagrangianFast(x,dx,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zero
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ba,...      % state to strain matrix
    Ktt,...     % geometric stiffness
    Mtt,...     % geometric inertia
    Zeta)        

% compute total length
n    = numel(x);
s    = 0;

Z1 = zeros(6,6+2*(n-1));
Z2 = zeros(n,3*n+1);
Z1(1:3,1:3) = Phi0;
Z1(1:3,4)   = p0;

%NLStiff = false; 

for ii = 1:(size(Th,3)/2)
    
    % first EL-diff eval
    [K1Z1,K1Z2] = LagrangianODEX(x, dx, Z1,...
        Th(:,:,2*ii-1), xia0(:,1,2*ii-1), Ba, Mtt, Ktt);
    
    % second EL-diff eval
    [K2Z1,K2Z2] = LagrangianODEX(x, dx, Z1 + (2/3)*ds*K1Z1,...
        Th(:,:,2*ii), xia0(:,1,2*ii), Ba, Mtt, Ktt);
    
    % update integrands
    s  = s  + ds;
    Z1 = Z1 + 0.25*ds*(K1Z1 + 3*K2Z1);
    Z2 = Z2 + 0.25*ds*(K1Z2 + 3*K2Z2);

end

% recover the kinematics entities
p   = Z1(1:3,4);
Phi = Z1(1:3,1:3);
B1  = Z1(1:6,5:5+n-1);
J   = Admapinv(Phi,p)*B1;

% recover the dynamics entities
M  = Z2(1:n,1:n);
C  = Z2(1:n,n+1:2*n);
K  = Z2(1:n,2*n+1:3*n);
G  = Z2(1:n,3*n+1);

Vg  = Z1(5,4);
Kin = Z1(6,4);

R = Zeta*K;

end

function [dZ1,dZ2] = LagrangianODEX(x,dx,Z1,...
    Theta,xia0,Ba,Mtt,Ktt)

n     = numel(x);
p_    = Z1(1:3,4);
Phi_  = Z1(1:3,1:3);
J_    = Z1(1:6,5:5+n-1);
Jt_   = Z1(1:6,6+n-1:6+2*(n-1));

%Theta_ = ThetaEval;%Model.ShpFnc(s);
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
dG = (Jg).'*(Ai*Mtt*[0;0;0;0;0;9.81e3]);

% compute (nonlinear stiffness)
dK = (BTh).'*Ktt*(BTh);

% compute grav. potential energy
dKe = 0.5*V.'*Mtt*V;
dVg = Mtt(4,4)*p_.'*[0;0;9.81e3];

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

end
%--------------------------------------------------------------------------
function y = isomSO3(x)
    x1 = x(1); x2 = x(2); x3 = x(3);
    y = [0, -x3, x2; x3, 0, -x1; -x2, x1, 0];
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