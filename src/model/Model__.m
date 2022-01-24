classdef Model
    
    properties (Access = public)
        Nlink;
        Sstep, Tstep;
        t, Tsim;
        q, q0;
        dq, dq0;
        lam, lam0;
        p, Phi, p0, Phi0;
        M, C, G, K, R, Vg, J;
        m0, l0, r0;
        Pi, Pihat, Y;
        hamil;
        Shp;
        tau, tau_;
        updatelaw;
        Payload;
    end
    
    properties (Access = private)
        N; S;
        Ktrue;
        ke, kb, kp, de, db;
        gam, mue, mub;
        dTaudq, dTauddq;
        ODE;
        MaxItr;
        Conv;
        Creep;
        Adaptive;
        Linewidth, Markersize;
    end
    
%--------------------------------------------------------------------------
methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------- Model Class
function obj = Model(N,varargin) 
    
    obj.Nlink  = N;
    obj.Sstep  = 40;
    obj.Tstep  = 1/50;
    obj.Tsim   = 10;
    obj.MaxItr = 3;
    obj.Conv   = 0.01;
    
    obj.ke  = [223.4, 174.0,-45.55];          % preset elongation stiff.   
    obj.kb  = [0.42292, 0.39552, -0.21293];   % preset bending stiff.
    
    obj.kp  = 0.75;                           % preset asymmetric stiff.
    obj.de  = 0.01;                           % preset elongation damp. 
    obj.db  = 1e-5;                           % preset bending damp.

    obj.mue = 0.523;
    obj.mub = 1.53e-2;
    obj.gam = [3.21e2, 5.22e-1, 26.356,...    % preset creeping para.
               0.00018193, 26.356, 0.00018193];

    obj.Phi0 = eye(3);
    obj.p0 = zeros(3,1);
    obj.m0 = ones(N,1)*0.05;
    obj.l0 = ones(N,1)*0.065;
    obj.r0 = ones(N,1)*0.05;
    
    obj.q0   = zeros(3*N,1);
    obj.dq0  = zeros(3*N,1);
    obj.lam0 = zeros(3*N,1);
    obj.tau  = @(mdl) zeros(3*N,1);
    
    obj.updatelaw = @(mdl) 0;
    
    obj.Shp      = chooseShapeFncfromLibary(obj,N);
    obj.ODE      = '.';
    obj.Creep    = false;
    obj.Adaptive = false;
    
    obj.Linewidth  = 2;
    obj.Markersize = 18;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    if ~obj.Creep, obj.gam = obj.gam*0; obj.lam0 = obj.lam0*0; end
    
    % build parameter vector
    obj = rebuildParameters(obj);
    obj.Pihat = obj.Pi;

end
%---------------------------------------------------------------------- get     
function varargout = get(Model,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Model.(varargin{ii});
        end
    else
        varargout = Model.(varargin);
    end
end   
%---------------------------------------------------------------------- set
function Model = set(Model,varargin)
    
    for ii = 1:2:length(varargin)
        if strcmp(varargin{ii},'r0') || strcmp(varargin{ii},'m0') || ...
                strcmp(varargin{ii},'l0') 
           if numel(varargin{ii+1}) == 1
                Model.(varargin{ii}) = varargin{ii+1}*ones(Model.Nlink,1);
           else
                Model.(varargin{ii}) = varargin{ii+1};
           end
        else 
            Model.(varargin{ii}) = varargin{ii+1};
        end
    end
    
    if Model.Creep
        Model.gam = [3.21e2, 5.22e-1, 26.356,...
                     0.0018193, 26.356, 0.0018193];
    end
    
    if ~isa(Model.ke,'function_handle') || ~isa(Model.kb,'function_handle')
        Model = rebuildParameters(Model);
    end
end
%------------------------------------------ set number of discrete elements 
function Model = setElements(Model,N), Model.Sstep = N; end
%-------------------------------------------- set the implicit solver freq.
function Model = setFrequency(Model,f), Model.Tstep = 1/f; 
Model.t = stepspace(0,Model.Tsim,Model.Tstep).';
end
%-------------------------------------------------------- set mass of links
function Model = setMass(Model,a) 
    if numel(a) == 1, Model.m0 = ones(Model.Nlink,1)*a; 
    else, Model.m0 = a(:);
    end
end
%------------------------------------------------------ set length of links
function Model = setLength(Model,a)
    if numel(a) == 1, Model.l0 = ones(Model.Nlink,1)*a; 
    else, Model.l0 = a(:);
    end
end
%------------------------------------------------------ set length of links
function Model = setRadius(Model,a)
    if numel(a) == 1, Model.r0 = ones(Model.Nlink,1)*a; 
    else, Model.r0 = a(:);
    end
end
%------------------------------------------------------ set length of links
function Model = setDamping(Model,a)
    if numel(a) == 1, Model.db = a; 
    else, Model.de = a(1); Model.db = a(2); 
    end
end
%------------------------------------------------------ set shape-functions
function Model = setShapeFunction(Model,fnc), Model.Shp = fnc; end
%------------------------------------------------------ set shape-functions
function Model = setLoad(Model,m)
    Model.Payload = m; 
    Model = rebuildParameters(Model);
    Model.Pihat = Model.Pi;
end
%------------------------------------------------------ set MB-controllers
function Model = setControl(Model,fnc)
    Model.tau = fnc; 
    [Khat,Dhat] = computeControlJacobians(Model);
    
    Model.dTaudq = Khat;
    Model.dTauddq = Dhat;
end
%---------------------------------------------------- simulates dyn. system
function [Model,treal] = simulate(Model)
    
Model.t  = stepspace(0,Model.Tsim,Model.Tstep);
Model.q  = []; Model.dq = [];

if isempty(Model.dTaudq)
    [Model.dTaudq,Model.dTauddq] = computeControlJacobians(Model);
end

if strcmp(Model.ODE,'.')
    [T, X, U, Pih, Hm, treal] = simulateSoftRobot(Model,...
        [Model.q0(:);Model.dq0(:);Model.lam0(:)]);
else
    [T, X, U, Pih, Hm, treal] = simulateSoftRobotODEMATLAB(Model,...
        [Model.q0(:);Model.dq0(:);Model.lam0(:)]);
end

% extracting data
Model.t     = T(:);
Model.q     = X(:,1:(3*Model.Nlink),:);
Model.dq    = X(:,3*Model.Nlink+1:6*Model.Nlink);  
Model.lam   = X(:,6*Model.Nlink:end);  
Model.tau   = U;
Model.hamil = Hm;
Model.Pihat = Pih;

end
%--------------------------------------------- compute end-effector pos/vel 
function [P, V] = computeEndEffector(Model,Q,dQ)
[P, ~, J_] = computeForwardKinematics(Model,Q);
P = P(end,:).';
V = J_*dQ(:);
end
%------------------------------- draws soft robot for current configuration
function [P] = show(Model,Q,col)
    
    [P, Np] = computeForwardKinematics(Model,Q);
    
    if nargin < 2
        col = 'b';
    end
    
    plt1 = plot3(P(:,1),P(:,2),P(:,3),'-','linewidth',...
        Model.Linewidth,'Color',col); hold on;
    plt2 = plot3(P(Np,1),P(Np,2),P(Np,3),'.',...
        'markersize',Model.Markersize,'Color',col);
    
    plt = {plt1,plt2};
end

end
%--------------------------------------------------------------------------
methods (Access = private) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------- implicit time-integration solver
function [Ts, X, U, PI, Hm, tt] = simulateSoftRobot(Model,z0)

Ts  = Model.t(:);
h   = mean(diff(Ts));
nd  = Model.Nlink*3;

z  = z0(:); 
X  = zeros(length(Ts)-1,nd*3);
U  = zeros(length(Ts)-1,nd);
Hm = zeros(length(Ts)-1,1);

X(1,:)  = z0(:).';
U(1,:)  = zeros(1,nd);
PI(1,:) = Model.Pihat(:).';

% pre-compute the controller jacobians for implicit solver
if isempty(Model.dTaudq)
    [Khat,Dhat] = computeControlJacobians(Model);
    
    Model.dTaudq = Khat;
    Model.dTauddq = Dhat;
end

tic;
disp('----------------------------------');
fprintf('* Computing SR dynamics... ');
progress('start')

for ii = 1:length(Ts)-1
    
    % assign states to temp. variable
    dW  = 1e3;
    w_  = z;
    itr = 1;
    
    while (dW > Model.Conv && itr <= Model.MaxItr)
        
        % compute flow field
        if itr == 1
            [f1, Model, Mi] = flow(Model,z,Ts(ii));
            dR = -h*f1;
        else
            [f2, Model, Mi] = flow(Model,w_,Ts(ii)+h);
            dR = w_ - z - 0.5*h*(f1 + f2);
        end
        
        % compute hessian
        H = buildHessian(Model,Mi);
        
        % hessian update iteration
        dw = stateUpdate(Model,H,dR);
        
        % update local state solution
        w_  = w_  + dw;
        itr = itr + 1;
        
        % compute convergence 
        dW = norm(dw(1:3*Model.Nlink));
    end
    
    % update states
    z = w_;
    
    % update parameters
    if Model.Adaptive
        dPi = Model.updatelaw(Model);
        Model.Pihat = Model.Pihat + h*dPi;
    end
   
    % write output data
    X(ii+1,:)  = z(:).';
    PI(ii+1,:) = Model.Pihat(:).';
    U(ii+1,:)  = Model.tau_(:).';
    Hm(ii)     = 0.5*z(nd+1:2*nd).'*Model.M*z(nd+1:2*nd) + ... 
                 0.5*z(1:nd).'*Model.K*z(1:nd) + Model.Vg;
    
    % update progress bar
    inc = round(100*(ii/(length(Ts)))-1);
    progress(inc,100);
end

Hm = [Hm; Hm(end)];

progress(100,100);
progress('end');
disp('----------------------------------');

tt = toc;
disp(['* Number of elem.  = ', num2str(Model.Sstep,3)]);
disp(['* Computation time = ', num2str(tt,3), ' s']);
disp(['* Computation freq = ', num2str(round(1/h)), ' Hz']);
disp(['* Real-time ratio  = ', num2str((Ts(end))/(tt),4), '']);

disp('----------------------------------');

    function [f, Model, Minv] = flow(Model,Z,T)
        
        n  = numel(Model.m0);
        Q  = Z(1:3*n);
        dQ = Z(3*n+1:6*n);
        Lm = Z(6*n+1:end);
       
        % compute Lagrangian entries
        [M_,C_,K_,Kh_,R_,G_,...
         p_,Phi_,J_,Vg_,G1,G2,Klm,Y_] = computeLagrangian(Model,Q,dQ);
        
        % overwrite dynamics
        Model.q  = Q;
        Model.dq = dQ;
        Model.M  = M_;
        Model.C  = C_;
        Model.R  = R_;
        Model.G  = G_;
        Model.p  = p_;
        Model.Phi = Phi_;
        Model.J  = J_;
        Model.Vg = Vg_;
        Model.t  = T;
        
        if Model.Adaptive
            Model.K     = Kh_;
            Model.Ktrue = K_;
            Model.Y     = Y_;
        else
            Model.K = K_;
        end
        
        % evaluate control action
        Model.tau_ = Model.tau(Model);
        
        if ~isempty(Model.Payload)
            f = [0;0;0;0;0;-Model.Payload*9.81];
            deltaM = Model.J.'*(adjointSE3inv(Phi_,p_)*f);
            Model.tau_ = Model.tau_ + deltaM;
        end
        
        % evaluate creep strains
        Delta = Klm*Lm;
        
        % pre-compute Minverse
        Minv = M_\eye(numel(Q));
        
        % flow field
        f = [dQ; ...
             Minv*(Model.tau_ + Delta - C_*dQ - K_*Q - R_*dQ - G_); ...
             -G1*Lm - G2*dQ];
    end
    
end
%--------------------------------------- forwards integration of kinematics
function [pp, id, J] = computeForwardKinematics(Model,x)
    
% compute total length
ee   = x(1:3:numel(x));
ltot = sum(Model.l0(:).*(ee(:) + 1));

ds  = ltot/(Model.Sstep*5);
s   = 0.0;
p_   = Model.p0;
Phi_ = Model.Phi0;
J   = zeros(6,numel(x));

ss = [];
pp = [];

for ii = 1:(Model.Sstep*5)
    
   [K1p,K1Phi,K1J] = ForwardKinematicODE(Model,s,...
                        x, Phi_, p_);
    
   [K2p,K2Phi,K2J] = ForwardKinematicODE(Model, s + (2/3)*ds,...
                        x, Phi_ + (2/3)*ds*K1Phi, p_ + (2/3)*ds*K1p);
    
   s    = s    + ds; 
   p_   = p_   + 0.25*ds*(K1p + 3*K2p);
   Phi_ = Phi_ + 0.25*ds*(K1Phi + 3*K2Phi);
   J    = J    + 0.25*ds*(K1J + 3*K2J);
   
   ss = [ss; s];
   pp = [pp; p_(:).'];
   
end

% transform Jacobian to body frame
Ai = adjointSE3inv(Phi_,p_);
J  = Ai*J;

id = round(linspace(1,(Model.Sstep*5),Model.Nlink+1));

end
%---------------------------------------------- compute Lagrangian entities 
function [M,C,K,Khat,R,G,p,Phi,J,Vg,Kg1,Kg2,Klm,Y] = computeLagrangian(Model,x,dx)
ndof = numel(Model.m0);

% compute total length
n    = numel(x);
ee   = x(1:3:n);
ltot = sum(Model.l0(:).*(ee(:) + 1));

ds   = ltot/Model.Sstep;
p    = Model.p0;
Phi  = Model.Phi0;
K    = zeros(n);
Khat = zeros(n);
Y    = zeros(n,7);
s    = 0;

Pvec = zeros(Model.Nlink*3,1);
for ii = 1:length(Model.l0)
    Pvec(3*ii-2) = Model.m0(ii);
    Pvec(3*ii-1) = Model.l0(ii);
    Pvec(3*ii) = Model.r0(ii);
end

Z1 = zeros(6,4+2*n);
Z2 = zeros(n,2*n+1);
Z1(1:3,1:3) = Phi;
Z1(1:3,4)   = p;

for ii = 1:Model.Sstep
    
    % get lengths and eval PCC shapefunction
    l  = Model.l0(:).*(x(1:3:numel(x)) + 1);
    S_ = Model.Shp(s,l);
    
    parvec = S_*Pvec;
    
    [K1Z1,K1Z2] = LagrangianODEX(Model,s, x, dx, S_, parvec, Z1);
    
    [K2Z1,K2Z2] = LagrangianODEX(Model,s+(2/3)*ds, x, dx, S_, parvec,...
        Z1 + (2/3)*ds*K1Z1);
    
    % update integrands
    s  = s  + ds;
    Z1 = Z1 + 0.25*ds*(K1Z1 + 3*K2Z1);
    Z2 = Z2 + 0.25*ds*(K1Z2 + 3*K2Z2);
    
    % compute stiffness matrix    
    if Model.Adaptive
        [Ktt, Ktthat] = nonlinearStiffnessMat(Model,s,x);
        
        Ytt = buildRegressor(Model,s,x);
        
        Y = Y + (1/(ndof*ltot))*ds*(S_.'*Ytt);
        K = K + (1/(ndof*ltot))*ds*(S_.'*Ktt*S_);
        Khat = Khat + (1/(ndof*ltot))*ds*(S_.'*Ktthat*S_);
    else
        [Ktt, ~] = nonlinearStiffnessMat(Model,s,x);
        K    = K + (1/(ndof*ltot))*ds*(S_.'*Ktt*S_);
        Khat = K;
    end
end

R  = kron(eye(ndof),diag([Model.de,Model.db,Model.db]));

M  = Z2(1:n,1:n);
C  = Z2(1:n,n+1:2*n);
G  = Z2(1:n,2*n+1);
Vg = Z1(5,4);

p   = Z1(1:3,4);
Phi = Z1(1:3,1:3);
B1  = Z1(1:6,5:5+n-1);
J   = adjointSE3inv(Phi,p)*B1;

Kg1  = kron(eye(ndof),diag([Model.gam(1),Model.gam(3),Model.gam(5)]));
Kg2  = kron(eye(ndof),diag([Model.gam(2),Model.gam(4),Model.gam(6)]));
Klm  = kron(eye(ndof),diag([Model.mue,Model.mub,Model.mub]));

if ~isempty(Model.Payload)
    
    dMtt = Model.Payload*[zeros(3,3), zeros(3,3); zeros(3,3), eye(3,3)]; 
    
    B2 = Z1(1:6,5+n:end);
    dJ = adjointSE3inv(Phi,p)*B2;
    dM = J.'*dMtt*J;
    dC = J.'*dMtt*dJ;
    
    M = M + dM;
    C = C + dC;
    
    f = [0;0;0;0;0;-9.81];
    deltaMY = J.'*(adjointSE3inv(Phi,p)*f);
    
    Y = [Y,deltaMY];
    
end

end
%-------------------------------------------------- forwards kinematics ODE
function [dp,dPhi,dJ] = ForwardKinematicODE(Model,s,x,Phi_,p_)

% get lengths 
l = Model.l0(:).*(x(1:3:numel(x)) + 1);

% construct geometric vectors
Jstar = [0,0,-1;0,1,0;0,0,0;0,0,0;0,0,0;1,0,0];
XI    = Jstar*Model.Shp(s,l)*x(:) + [0,0,0,0,0,1].';

U     = XI(4:6);
Gamma = XI(1:3);

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*skew(Gamma);
A    = adjointSE3(Phi_,p_);
dJ   = A*Jstar*Model.Shp(s,l);

end
%---------------------------- Lagrangian Matrix-Differential Equation (MDE)
function [dZ1,dZ2] = LagrangianODEX(~,~,x,dx,S_,par_,Z1)
n     = numel(x);
m0_   = par_(1);
l0_   = par_(2);
r0_   = par_(3);

p_    = Z1(1:3,4);
Phi_  = Z1(1:3,1:3);
J_    = Z1(1:6,5:5+n-1);
Jt_   = Z1(1:6,6+n-1:6+2*(n-1));

% construct geometric vectors
Jstar = [0,0,-1;0,1,0;0,0,0;0,0,0;0,0,0;1,0,0];
XI    = Jstar*S_*x;

Gamma = XI(1:3);
U     = XI(4:6) + [0;0;1];

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*skew(Gamma);

A   = adjointSE3(Phi_,p_);
Ai  = adjointSE3inv(Phi_,p_);

% build jacobian
Jg  = Ai*J_;
Jgt = Ai*Jt_;

V  = Jg*dx;
Vw = V(1:3);
Vs = V(4:6);

a   = adjointse3(Vw,Vs);

dJ  = A*Jstar*S_;
dJt = A*a*Jstar*S_;

% build dynamic matrices
l   = l0_*(U(3));

msl = m0_/l;

Ixx = (r0_^2)/(4*l); 
Iyy = (r0_^2)/(4*l); 
Izz = (r0_^2)/(2*l); 

Is  = diag([Ixx,Iyy,Izz]); 

Mtt = zeros(6,6);
Mtt(1:3,1:3) = msl*Is;
Mtt(4:6,4:6) = msl*eye(3);
%    
% compute inertia, coriolis, gravity
dM = (Jg).'*Mtt*Jg;
dC = (Jg).'*((Mtt*a - a.'*Mtt)*Jg  + Mtt*Jgt);
dG = (Jg).'*msl*(Ai*[0;0;0;0;0;9.81]);

% compute grav. potential energy
dVg = msl*p_.'*[0;0;9.81];

dZ1 = zeros(6,4+2*n);
dZ1(1:3,1:3)             = dPhi;
dZ1(1:3,4)               = dp;
dZ1(1:6,5:5+n-1)         = dJ;
dZ1(1:6,6+n-1:6+2*(n-1)) = dJt;
dZ1(5,4)                 = dVg;

dZ2 = zeros(n,2*n+1);
dZ2(1:n,1:n)     = dM;
dZ2(1:n,n+1:2*n) = dC;
dZ2(1:n,2*n+1)   = dG;

end
%----------------------------------------- compute nonlinear stiffness mat.
function Model = rebuildParameters(Model)
if isempty(Model.Payload)
Model.Pi = [Model.ke(1),Model.ke(2),Model.ke(3),...
    Model.kb(1),Model.kb(2),Model.kb(3),Model.kp].';
else
Model.Pi = [Model.ke(1),Model.ke(2),Model.ke(3),...
    Model.kb(1),Model.kb(2),Model.kb(3),Model.kp,Model.Payload].';    
end
end
%--------------------------------------- updated Hessian with dyn. residual
function dr = stateUpdate(Model, H, dR)
dr = -(-(1/2)*Model.Tstep*H + eye(size(H,1),size(H,1)))\dR;
end
%---------------------------------------------- compute Hessian approximate
function DF = buildHessian(Model,varargin)
    
ndof = numel(Model.m0);
G1   = kron(eye(ndof),diag([Model.gam(1),Model.gam(3),Model.gam(5)]));
G2   = kron(eye(ndof),diag([Model.gam(2),Model.gam(4),Model.gam(6)]));
Klm  = kron(eye(ndof),diag([Model.mue,Model.mub,Model.mub]));

Minv = varargin{1};
n    = size(Model.q,1);

if isempty(Model.Ktrue)
   KK = Model.K;
else
   KK = Model.Ktrue;
end

DF                      = zeros(3*n,3*n);
DF(1:n,n+1:2*n)         = eye(n);
DF(n+1:2*n,1:n)         = -Minv*(KK - Model.dTaudq);
DF(n+1:2*n,n+1:2*n)     = -Minv*(Model.R + Model.C - Model.dTauddq);
DF(n+1:2*n,2*n+1:3*n)   = -Minv*Klm;
DF(2*n+1:3*n,n+1:2*n)   = -G2;
DF(2*n+1:3*n,2*n+1:3*n) = -G1;

end
%----------------------------------------- compute nonlinear stiffness mat.
function [Ktt, Ktthat] = nonlinearStiffnessMat(Model,s,x)
l  = Model.l0(:).*(x(1:3:numel(x)) + 1);
S_ = Model.Shp(s,l);

Q    = S_*x;  
parL = S_*[Model.l0(:);Model.l0(:);Model.l0(:)]; 
eps  = parL(2)*Q(1);
beta = parL(2)*(sqrt(Q(2)^2 + Q(3)^2));
phi  = atan2(Q(3),Q(2));

faxi = @(x,a) 0.5*a*(sin(3*(x)) + 1) + 1;
Fax = faxi(phi,Model.Pi(7)*beta);

if ~isa(Model.ke,'function_handle')
    KEE = Model.Pi(1) + Model.Pi(2)*(tanh(Model.Pi(3)*eps)^2 - 1);
else
    KEE = Model.ke(Q);
end

if ~isa(Model.kb,'function_handle')
    KBB = Fax*(Model.Pi(4) + Model.Pi(5)*(tanh(Model.Pi(6)*beta)^2 - 1)); 
else
    KBB = Model.kb(Q);
end

StiffnessCorrection = (numel(x)/3);

Ktt = StiffnessCorrection*...
    diag([KEE,KBB*parL(2),KBB*parL(2)]);       

Fax = faxi(phi,Model.Pihat(7)*beta);
KEE_ = Model.Pihat(1) + Model.Pihat(2)*(tanh(Model.Pihat(3)*eps)^2 - 1);
KBB_ = Fax*(Model.Pihat(4) + Model.Pihat(5)*(tanh(Model.Pihat(6)*beta)^2 - 1));

Ktthat = StiffnessCorrection*...
    diag([KEE_,KBB_*parL(2),KBB_*parL(2)]);   

end
%------------------------------------------------ compute Regressor matrix
function Y = buildRegressor(Model,s,x)
l  = Model.l0(:).*(x(1:3:numel(x)) + 1);
S_ = Model.Shp(s,l);

parL = S_*[Model.l0(:);Model.l0(:);Model.l0(:)]; 
Q    = S_*x;  
eps  = parL(2)*Q(1);
beta = parL(2)*(sqrt(Q(2)^2 + Q(3)^2));
phi  = atan2(Q(3),Q(2));

w    = 1/(2*pi/6);
faxi = @(x,a)-(0.5*cos(w*pi*(x) + pi/2)+0.5)*a+(1+a);    

L     = parL(2)*(1+Q(1));
Faxi  = faxi(phi,Model.Pihat(7)*beta);
alpha = (numel(x)/3);

Y1 = alpha;
Y2 = alpha*(tanh(Model.Pihat(3)*eps)^2 - 1);

Y4 = Faxi*alpha*L;
Y5 = Faxi*alpha*L*(tanh(Model.Pihat(6)*beta)^2 - 1);

YY1  = [Y1;0;0]*Q(1);
YY2  = [Y2;0;0]*Q(1);
YY3  = [0;0;0];
YY4  = [0;Y4*Q(2);Y4*Q(3)];
YY5  = [0;Y5*Q(2);Y5*Q(3)];
YY6  = [0;0;0];
YY7  = [0;0;0];

Y = [YY1,YY2,YY3,YY4,YY5,YY6,YY7];
  
end
%------------------------------------- select PCC shapefunction from libary
function shp = chooseShapeFncfromLibary(~,n)
   if n == 1, shp = @(x,l) PieceWiseFunction1(x,l);
   elseif n == 2, shp = @(x,l) PieceWiseFunction2(x,l);
   elseif n == 3, shp = @(x,l) PieceWiseFunction3(x,l);
   elseif n == 4, shp = @(x,l) PieceWiseFunction4(x,l);
   elseif n == 5, shp = @(x,l) PieceWiseFunction5(x,l);    
   elseif n == 6, shp = @(x,l) PieceWiseFunction6(x,l);
   elseif n == 7, shp = @(x,l) PieceWiseFunction7(x,l);
   elseif n == 8, shp = @(x,l) PieceWiseFunction8(x,l);     
   elseif n == 9, shp = @(x,l) PieceWiseFunction9(x,l);     
   else, error(['PCC shape function order available until n=9. \n If more links'...
           ,' are needed, please modify the PieceWiseFunction_.m files under src/pwf']);
   end
end
%----------------------------------- (pre)-computes the controller jacobian
function [Kt, Dt] = computeControlJacobians(Model)
n   = length(Model.q0);
Q0  = Model.q0(:);    
dQ0 = Model.dq0(:);    
    
% compute Lagrangian entries
[M_,C_,K_,Kh_,R_,G_,p_,Phi_,J_] = computeLagrangian(Model,Q0,dQ0);

% overwrite dynamics
Model.q  = Q0; Model.dq = dQ0;
Model.M  = M_; Model.C  = C_; Model.R  = R_; Model.G  = G_; 
Model.p  = p_; Model.Phi = Phi_; Model.J  = J_;
Model.t  = 0;

if Model.Adaptive
    Model.K = Kh_;
else
    Model.K = K_;
end
    
epsilon = 1e-3;
delta   = epsilon*eye(n);

Tau0 = Model.tau(Model);

Kt = zeros(n);
Dt = zeros(n);

for ii = 1:n
    Model.q = Q0 + delta(:,ii);
    Tau_ = Model.tau(Model);
    Kt(:,ii) = (Tau_ - Tau0)/epsilon;
end

Model.q  = Q0; 

for ii = 1:n
    Model.dq = dQ0 + delta(:,ii);
    Tau_ = Model.tau(Model);
    Dt(:,ii) = (Tau_ - Tau0)/epsilon;
end

end

end
end
