classdef Model
    
    properties (Access = public)
        Shapes;
        Inertia;
        Tsim;
        Tstep;

        q, dq, t;
        q0, dq0, Phi0, p0;
        Log;
   
        tau, tau_;
    end
    
    properties (Access = private)
        N; S;
        dTaudq, dTauddq;
        Theta;
        Xi0;
        
        MexSolver = true;
        
        ResidualNorm;
        MaxIteration;
        Conv;
        Adaptive;
        
        Linewidth;
        Markersize;
        
        Ba;
        ShpFnc;
    end
    
%--------------------------------------------------------------------------
methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------- Model Class
function obj = Model(Shapes,varargin) 
    
    obj.Shapes = Shapes;
    obj.Tstep  = 1/60;
    obj.Tsim   = 10;
    obj.Ba     = Shapes.Ba;
    obj.ShpFnc = Shapes.Theta;
    
    obj.MaxIteration = 3;
    obj.ResidualNorm = 0.01;

    G0 = Shapes.get('g0');
    obj.Phi0 = G0(1:3,1:3); 
    obj.p0   = G0(5:7).';    
    obj.q0   = zeros(Shapes.NDim,1) + 1e-3*rand(Shapes.NDim,1);
    obj.dq0  = zeros(Shapes.NDim,1);

    obj.tau  = @(mdl) zeros(Shapes.NDim,1);
    
    obj.Linewidth  = 4;
    obj.Markersize = 25;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end

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
        if strcmp(varargin{ii},'r0')

        else 
            Model.(varargin{ii}) = varargin{ii+1};
        end
    end
    
end
%---------------------------------------------------- simulates dyn. system
function Model = simulate(Model)
    
Model.t   = 0:Model.Tstep:Model.Tsim;
Model.q   = []; 
Model.dq  = [];

Model.Log    = struct;
Model.Log.q  = Model.q0;
Model.Log.dq = 0.1*Model.q0;
Model.Log.EL = struct;

% get evals
Model.Theta = Model.Shapes.get('ThetaEval');
Model.Xi0 = Model.Shapes.get('Xi0Eval');


if isempty(Model.dTaudq)
    Model.Log.t = 0;
    [Model.dTaudq,Model.dTauddq] = computeControlJacobians(Model);
end

[T, X, U, Km, Ue, Ug] = simulateSoftRobot(Model,...
    [Model.q0(:); Model.dq0(:)]);

% extracting data
Model.Log.t   = T(:);
Model.Log.q   = X(:,1:Model.Shapes.NDim,:);
Model.Log.dq  = X(:,2*Model.Shapes.NDim+1:Model.Shapes.NDim);  
Model.Log.tau = U;

Model.Log.Kin = Km;
Model.Log.Vg  = Ug;
Model.Log.Psi = Ue;

end
%---------------------------------------------------- simulates dyn. system
function Model = computeEL(Model,Q,varargin)

if isempty(varargin), dQ = Q*0; end    
    
% compute Lagrangian entries
[M_,C_,K_,R_,G_,...
p_,Phi_,J_,Vg_,Kin_] = computeLagrangian(Model,Q,dQ);
    
% overwrite dynamics
%Model.Log.t  = T;
Model.Log.q  = Q;
Model.Log.dq = dQ;
Model.Log.p   = p_;
Model.Log.Phi = Phi_;

Model.Log.EL.M = M_;
Model.Log.EL.C = C_;
Model.Log.EL.R = diag(diag(R_));
Model.Log.EL.K = diag(diag(K_));
Model.Log.EL.G = G_;
Model.Log.EL.J = J_;

Model.Log.Vg  = Vg_;
Model.Log.Kin = Kin_;

end
%------------------------------------------------- plot curve configuration
function [P] = show(Model,Q,col,varargin)
    
    % solve forward kinematics on [0,L0]
    %[P, Np] = computeForwardKinematics(Model,Q(:));
    
    P = Model.Shapes.string(Q(:));
    
    if nargin < 2
        col = col(1);
    end
    
    if ~isempty(varargin)
       st = '--'; 
    else
        st = '-';
    end
       
    % plot spatial curve
    plot(P(:,1),P(:,3),'-','linewidth',...
         Model.Linewidth,'Color',col,'LineStyle',st); hold on;
    % plot interconnection links 
    plot(P([1,end],1),P([1,end],3),'.',...
         'markersize',Model.Markersize,'Color',col);
end
%--------------------------------------------- compute end-effector pos/vel 
function [p, eta] = endeffector(Model,Q,dQ)
[p, ~, Jacob] = computeForwardKinematics(Model,Q(:));
p   = p(end,:).';
eta = Jacob*dQ(:);
end

end
%--------------------------------------------------------------------------
methods (Access = private) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------- implicit time-integration solver
function [Ts, X, U, Kin, Ue, Ug] = simulateSoftRobot(Model,z0)

Ts  = Model.t(:);
h   = mean(diff(Ts));
nd  = Model.Shapes.NDim;

z   = z0(:); 
X   = zeros(length(Ts),nd*2);
U   = zeros(length(Ts),nd);
Kin = zeros(length(Ts),1);
Ue  = zeros(length(Ts),1);
Ug  = zeros(length(Ts),1);

X(1,:)  = z0(:).';
U(1,:)  = zeros(1,nd);

tic;
disp('----------------------------------');
fprintf('* Computing SR dynamics... ');
progress('start')

for ii = 1:length(Ts)-1
    
    % assign states to temp. variable
    dW  = 1e3;
    w_  = z;
    itr = 1;
    
    while (dW > Model.ResidualNorm && itr <= Model.MaxIteration)
        
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
        dW = norm(dw(1:Model.Shapes.NDim));
    end
    
    % update states
    z = w_;
       
    % write output data
    X(ii+1,:) = z(:).';
    U(ii+1,:) = Model.tau_(:).';
    Kin(ii+1) = 0.5*z(nd+1:2*nd).'*Model.Log.EL.M*z(nd+1:2*nd);
    Ue(ii+1)  = 0.5*z(1:nd).'*Model.Log.EL.K*z(1:nd);
    Ug(ii+1)  = Model.Log.Vg;

    % update progress bar
    inc = round(100*(ii/(length(Ts)))-1);
    progress(inc,100);
end

progress(100,100);
progress('end');
disp('----------------------------------');

tt = toc;
disp(['* Number of elem.  = ', num2str(Model.Shapes.NNode,3)]);
disp(['* Computation time = ', num2str(tt,3), ' s']);
disp(['* Computation freq = ', num2str(round(1/h)), ' Hz']);
disp(['* Real-time ratio  = ', num2str((Ts(end))/(tt),4), '']);

disp('----------------------------------');

    function [f, Model, Minv] = flow(Model,Z,T)
        
        n  = Model.Shapes.NDim;
        Q  = Z(1:n);
        dQ = Z(n+1:2*n);
       
        % compute Lagrangian entries
        if ~Model.MexSolver
            [M_,C_,K_,R_,G_,...
                p_,Phi_,J_,Vg_,Kin_] = computeLagrangian(Model,Q,dQ);
        else
            [M_,C_,K_,R_,G_,...
                p_,Phi_,J_,Vg_,Kin_] = computeLagrangianFast_mex(...
                Q,dQ,... 
                Model.Shapes.ds,...   
                Model.p0,... 
                Model.Phi0,...
                Model.Xi0,... 
                Model.Theta,...
                Model.Shapes.Ba,... 
                Model.Shapes.Ktt,...
                Model.Shapes.Mtt,...     
                Model.Shapes.Zeta);      
        end
            
        % overwrite dynamics
        Model.Log.t  = T;
        Model.Log.q  = Q;
        Model.Log.dq = dQ;
        Model.Log.p   = p_;
        Model.Log.Phi = Phi_;
        
        Model.Log.EL.M = M_;
        Model.Log.EL.C = C_;
        Model.Log.EL.R = R_;
        Model.Log.EL.K = K_;
        Model.Log.EL.G = G_;
        Model.Log.EL.J = J_;

        Model.Log.Vg  = Vg_;
        Model.Log.Kin = Kin_;
        
        % evaluate control action
        Model.tau_ = Model.tau(Model);
        
        % pre-compute Minverse
        Minv = M_\eye(numel(Q));
        
        % flow field
        f = [dQ; ...
             Minv*(Model.tau_ - C_*dQ - K_*Q - R_*dQ - G_)];
    end
    
end
%--------------------------------------- forwards integration of kinematics
function [pp, id, J] = computeForwardKinematics(Model,x)
    
% compute total length
ds  = Model.Shapes.ds;
s    = 0.0;
p_   = Model.p0;
Phi_ = Model.Phi0;
J    = zeros(6,numel(x));

pp = p_';

% numerical trapzoid solver over [0,L0]
for ii = 1:Model.Shapes.NNode-1
    
   s = Model.Shapes.Sigma(ii);
    
   % first trapzoidal step 
   [K1p,K1Phi,K1J] = ForwardKinematicODE(Model,s,...
                        x, Phi_, p_);
   % second trapzoid step
   [K2p,K2Phi,K2J] = ForwardKinematicODE(Model, s + (2/3)*ds,...
                        x, Phi_ + (2/3)*ds*K1Phi, p_ + (2/3)*ds*K1p); 
                    
   p_   = p_   + 0.25*ds*(K1p + 3*K2p);
   Phi_ = Phi_ + 0.25*ds*(K1Phi + 3*K2Phi);
   J    = J    + 0.25*ds*(K1J + 3*K2J);
   
   pp = [pp;p_.'];   
   
end

% transform Jacobian to body frame
J  = adjointSE3inv(Phi_,p_)*J;
id = round(linspace(1,Model.Shapes.NNode,2));

end
%---------------------------------------------- compute Lagrangian entities 
function [M,C,K,R,G,p,Phi,J,Vg,Kin] = computeLagrangian(Model,x,dx)

% compute total length
n    = numel(x);
ds   = Model.Shapes.ds;
p    = Model.p0;
Phi  = Model.Phi0;
s    = 0;

% pre-computed Theta and Xi0;
Th = Model.Theta;
Xi = Model.Xi0;

Z1 = zeros(6,6+2*(n-1));
Z2 = zeros(n,3*n+1);
Z1(1:3,1:3) = Phi;
Z1(1:3,4)   = p;

if isa(Model.Shapes.Ktt,'function_handle')
   NLStiff = true; 
else
   NLStiff = false; 
end

for ii = 1:Model.Shapes.NNode
    
    % first EL-diff eval
    [K1Z1,K1Z2] = LagrangianODEX(Model,Th(:,:,2*ii-1),Xi(:,1,2*ii-1),...
        x, dx, Z1,NLStiff);
    
    % second EL-diff eval
    [K2Z1,K2Z2] = LagrangianODEX(Model,Th(:,:,2*ii), Xi(:,1,2*ii),...
        x, dx, Z1 + (2/3)*ds*K1Z1, NLStiff);
    
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

R = Model.Shapes.Zeta*K;

end
%-------------------------------------------------- forwards kinematics ODE
function [dp,dPhi,dJ] = ForwardKinematicODE(Model,s,x,Phi_,p_)

% construct geometric vectors
Theta = Model.Shapes.Theta(s);
XI    = Model.Shapes.Ba*Theta*x + Model.Shapes.xia0;

U     = XI(4:6);
Gamma = XI(1:3);

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*skew(Gamma);
A    = adjointSE3(Phi_,p_);
dJ   = A*Model.Shapes.Ba*Theta;
end

%---------------------------- Lagrangian Matrix-Differential Equation (MDE)
function [dZ1,dZ2] = LagrangianODEX(Model,Theta_,Xi0_,x,dx,Z1,NLStiff)

n     = numel(x);
p_    = Z1(1:3,4);
Phi_  = Z1(1:3,1:3);
J_    = Z1(1:6,5:5+n-1);
Jt_   = Z1(1:6,6+n-1:6+2*(n-1));

%Theta_ = ThetaEval;%Model.ShpFnc(s);
XI = Model.Ba*Theta_*x + Xi0_;

Gamma = XI(1:3);
U     = XI(4:6);

% build forward kin - position
dp   = Phi_*U;
dPhi = Phi_*isomSO3(Gamma);

A   = Admap(Phi_,p_);
Ai  = Admapinv(Phi_,p_);

% build jacobian
Jg  = Ai*J_;
Jgt = Ai*Jt_;

V    = Jg*dx;
adXi = admap(XI);
adV  = admap(V);

BTh = Model.Shapes.Ba*Theta_;

dJ  = A*BTh;
dJt = A*adV*BTh;
Mtt = Model.Shapes.Mtt;

if ~NLStiff
    Ktt = Model.Shapes.Ktt;
else
    Ktt = Model.Shapes.Ktt(XI);
end

%deta = -adXi*V + BTh*dx;

% compute inertia, coriolis, gravity
dM = (Jg).'*Mtt*Jg;
dC = (Jg).'*((Mtt*adV - adV.'*Mtt)*Jg  + Mtt*Jgt);
dG = (Jg).'*(Ai*Mtt*[0;0;0;0;0;-9.81e3]);

% compute (nonlinear stiffness)
%Ktt = nonlinearStiffnessMat(Model,s,x);
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
dZ2(1:n,1:n)       = dM;
dZ2(1:n,n+1:2*n)   = dC;
dZ2(1:n,2*n+1:3*n) = dK;
dZ2(1:n,3*n+1)     = dG;

end
%--------------------------------------- updated Hessian with dyn. residual
function dr = stateUpdate(Model, H, dR)
dr = -(-(1/2)*Model.Tstep*H + eye(size(H,1),size(H,1)))\dR;
end
%---------------------------------------------- compute Hessian approximate
function DF = buildHessian(Model,varargin)

Minv = varargin{1};
n    = size(Model.Log.q,1);

DF                      = zeros(2*n,2*n);
DF(1:n,n+1:2*n)         = eye(n);
DF(n+1:2*n,1:n)         = -Minv*(Model.Log.EL.K - Model.dTaudq);
DF(n+1:2*n,n+1:2*n)     = -Minv*(Model.Log.EL.R + Model.Log.EL.C - Model.dTauddq);

end
%--------------------------------------------------------- show solver info
function showInformation(Model)   

fprintf('--------------------------------------------------------------\n');  
fprintf('* Element = %i \n',NodeNum);
fprintf('* Max iteration = %i \n', Fem.MaxIteration);
fprintf('* Solver time horizon = %i \n', Nodel.TimeEnd);
fprintf('* Solver time step    = %i1.1e \n', Model.TimeStep);
showMaterialInfo(Fem)
fprintf('--------------------------------------------------------------\n');

end
%----------------------------------------- compute nonlinear stiffness mat.
function Ktt = nonlinearStiffnessMat(Model,s,x)

Jt = Model.Inertia.Jtt;

E0 = (Model.E);
G0 = (Model.E)/(2*(1+Model.Nu));
  
QE = diag([E0,G0,G0]);
QR = diag([G0,E0,E0]);

Theta1 = Model.Shapes.Phi(s);
Ba    = Model.Shapes.Ba;
xi    = Ba*Theta1*x;

%diag([G0*J11,E0*J22,E0*J33,E0*A,G0*A,G0*A]);
Ktt = blkdiag(QR*Jt,QE*eye(3));
Kappa = norm(abs(xi(1:3)));
Gamma = norm(abs(xi(4:6)));

Ktt = Ktt*(1.0 + 1.35*(tanh(-16200*Kappa)^2));

end
%----------------------------------- (pre)-computes the controller jacobian
function [Kt, Dt] = computeControlJacobians(Model)
n   = length(Model.q0);
Q0  = Model.q0(:);    
dQ0 = Model.dq0(:);    
    
% compute Lagrangian entries
[M_,C_,K_,R_,G_,p_,Phi_,J_] = computeLagrangian(Model,Q0,dQ0);

% overwrite dynamics
Model.q  = Q0; 
Model.dq = dQ0;
Model.Log.EL.M  = M_; 
Model.Log.EL.C  = C_; 
Model.Log.EL.R  = R_; 
Model.Log.EL.G  = G_; 
Model.Log.EL.K  = K_; 
Model.Log.EL.J  = J_;

Model.Log.p   = p_; 
Model.Log.Phi = Phi_; 

Model.t  = 0;

epsilon = 1e-3;
delta   = epsilon*eye(n);

Tau0 = Model.tau(Model);

Kt = zeros(n);
Dt = zeros(n);

% finite difference for tau(q(t),.)
for ii = 1:n
    Model.q = Q0 + delta(:,ii);
    Tau_ = Model.tau(Model);
    Kt(:,ii) = (Tau_ - Tau0)/epsilon;
end

Model.q  = Q0; 

% finite difference for tau(q(t),.)
for ii = 1:n
    Model.dq = dQ0 + delta(:,ii);
    Tau_ = Model.tau(Model);
    Dt(:,ii) = (Tau_ - Tau0)/epsilon;
end

end

end
end
