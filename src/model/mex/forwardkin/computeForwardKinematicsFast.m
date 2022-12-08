%#codegen
function [gtmp,Jtmp,vtmp] = computeForwardKinematicsFast(x,dx,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zero
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ba)         % state to strain matrix        

% compute total length
Ns   = (size(Th,3)/2);
n    = numel(x);
s    = 0;

Z1 = zeros(6,5+n-1);
Z1(1:3,1:3) = Phi0;
Z1(1:3,4)   = p0;

gtmp = zeros(4,4,Ns); 
Jtmp = zeros(6,n,Ns); 
vtmp = zeros(6,1,Ns); 

for ii = 1:Ns
    
    % first EL-diff eval
    [K1Z1] = ForwardODEX(x, Z1,...
        Th(:,:,2*ii-1), xia0(:,1,2*ii-1), Ba);
    
    % second EL-diff eval
    [K2Z1] = ForwardODEX(x, Z1 + (2/3)*ds*K1Z1,...
        Th(:,:,2*ii), xia0(:,1,2*ii), Ba);
    
    % update integrands
    s  = s  + ds;
    Z1 = Z1 + 0.25*ds*(K1Z1 + 3*K2Z1);
    
    % recover the kinematics entities
    p   = Z1(1:3,4);
    Phi = Z1(1:3,1:3);
    B1  = Z1(1:6,5:5+n-1);
    
    gtmp(:,:,ii) = SE3(Phi,p);
    Jtmp(:,:,ii) = Admap(Phi,[0;0;0])*Admapinv(Phi,p)*B1;
    vtmp(:,:,ii) = Jtmp(:,:,ii)*dx;
end

end

%--------------------------------------------------------------------------
function [dZ1] = ForwardODEX(x,Z1,Theta,xia0,Ba)

n     = numel(x);
p_    = Z1(1:3,4);
Phi_  = Z1(1:3,1:3);

BTh = Ba*Theta;
XI  = BTh*x + xia0;

Gamma = XI(1:3);
U     = XI(4:6);

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*isomSO3(Gamma);
A    = Admap(Phi_,p_);
dJ   = A*BTh;

dZ1              = zeros(6,5+n-1);
dZ1(1:3,1:3)     = dPhi;
dZ1(1:3,4)       = dp;
dZ1(1:6,5:5+n-1) = dJ;

end
%--------------------------------------------------------------------------
function Y = SE3(R,p)
    Y          = zeros(4);
    Y(4,4)     = 1;
    Y(1:3,1:3) = R;
    Y(1:3,4)   = p;
end
%--------------------------------------------------------------------------
function y = isomSO3(x)
    x1 = x(1); x2 = x(2); x3 = x(3);
    y = [0, -x3, x2; 
         x3, 0, -x1; 
         -x2, x1, 0];
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