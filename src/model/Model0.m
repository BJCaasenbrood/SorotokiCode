classdef Model

    properties (Access = public)
        Links;
        Dof;
        NDof;
        NLink = 1;
        tspan = 1;
        Table;
        g;
    end
    
    properties (Access = private)
        Nq;
        q; 
        dq;
        ddq;
        t;
        E = 250;
        nu = 0.495;
        mu = .1;
        xia0 = [0,0,0,1,0,0].';
        Phi;
        Ba; Bc; 
        Mee; Kee; Dee; 
        NModal = 5;
        Grav = 9.81;
        Radius = 1e-2;
        Density = 1;
        Length = 1;
        Length0;
        P1; P2; P3;
    end
    
%--------------------------------------------------------------------------
methods  
%-------------------------------------------------------------- Model Class
function obj = Model(Table,varargin) 
    
    obj.Table = Table;
    obj.NLink = size(Table,1);
    
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
        Model.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- set
function Model = GenerateModel(Model)
    Model.Length0 = Model.Length;
    Model = GenerateCosserat(Model);
end

%---------------------------------------------------------------------- set
function Model = solve(Model)
    
T = linspace(0,Model.tspan,1e2);
q0 = zeros(Model.NModal*Model.NDof,1);

%options = odeset('Events', @myEvent);
opt = odeset('RelTol',1e-3,'AbsTol',1e-5);
[ts,ys] = ode45(@(t,x) KinematicODE(Model,t,x),T,q0,opt);

Model.g = ys;
Model.t = ts;

end

%----------------------------------------------------------------- simulate
function Model = simulate(Model)
    
dt = 1e-2;
T = [0 Model.tspan];
q0 = zeros(Model.NModal*Model.NDof,1);
dq0 = zeros(Model.NModal*Model.NDof,1);

%options = odeset('Events', @myEvent);
opt = odeset('RelTol',1e-3,'AbsTol',1e-3);
[ts,ys] = ode23tb(@(t,x) DynamicODE(Model,t,x),T,[q0 dq0],opt);

Model.g = ys;
Model.t = ts;

end

%--------------------------------------------------------------------- show
function h = show(Model)
   
N = Model.Nq;
X = linspace(0,1,200);

Model.q(:,1) = Model.g(end,1:N).';
Model.dq(:,1) = 0*Model.q(:,1);
Model.ddq(:,1) = 0*Model.q(:,1);
ee = Model.Ba*Model.Phi(Model.Length)*Model.q;
Model.Length = ee(4);
g0 = [1,0,0,0,0,0,0];
eta0 = [0,0,0,0,0,0];
deta0 = [0,0,0,0,0,0];
[~,yf] = ode45(@(t,x) ForwardODE(Model,t,x),X,[g0 eta0 deta0]);
G = yf(:,1:7);

figure(101);
plot3(G(:,5),G(:,6),G(:,7),'linewidth',2,'Color',col(1));
axis equal;
    
end

%------------------------------------------------------------ show 3D model
function h = showModel(Model)
    
figure(101); hold all;
N = Model.Nq;
Model.Length = Model.Length0;
X = linspace(0,1,200);

msh = Gmodel('SoftActuatorRedux.stl');
mshgr = Gmodel('SoftGripperRedux.stl');

%% set texture
msh.Texture = grey;
mshgr.Texture = grey;

msh = msh.bake();
mshgr = mshgr.bake();
msh = msh.render();
mshgr = mshgr.render();

LinkID = knnsearch(X(:),msh.Node(:,3));
axis equal;
drawnow;
axis([-.5 .5 -.5 .5 -1 0])
view(40,5);
msh.update();
mshgr.update();

FPS = round((1/12)/(mean(diff(Model.t))));

for ii = 1:FPS:length(Model.t)

    msh.resetNode();
    mshgr.resetNode();

    Model.q(:,1) = Model.g(ii,1:N).';
    Model.dq(:,1) = 0*Model.q(:,1);
    Model.ddq(:,1) = 0*Model.q(:,1);
    ee = Model.Ba*Model.Phi(Model.Length)*Model.q;
    Model.Length = ee(4);
    g0 = [1,0,0,0,0,0,0];
    eta0 = [0,0,0,0,0,0];
    deta0 = [0,0,0,0,0,0];
    [~,yf] = ode45(@(t,x) ForwardODE(Model,t,x),X,[g0 eta0 deta0]);

    SweepSE3 = yf(:,1:7);
    
    msh = Blender(msh,'Rotate',{'z',30});
    mshgr = Blender(mshgr,'Rotate',{'z',30});
    msh = Blender(msh,'Sweep', {LinkID,SweepSE3});
    msh = Blender(msh,'Rotate',{'x',-180});
    mshgr = Blender(mshgr,'SE3',yf(end,1:7));
    mshgr = Blender(mshgr,'Rotate',{'x',-180});
    
    mshgr.updateNode();
    msh.updateNode();
    msh.update();
    mshgr.update();
    axis(1.5*[-.5 .5 -.5 .5 -1 0]);
    drawnow;
end
    
end

end
%--------------------------------------------------------------------------
methods (Access = private)

%---------------------------------------------------------------------- set
function Model = GenerateCosserat(Model)
I6 = eye(6);
set = 1:6;
xa = set(logical(Model.Table));
xc = setdiff(set,xa);
Model.Ba = I6(:,xa);
Model.Bc = I6(:,xc);
Model.NDof = size(Model.Ba,2);

Model.Phi = ShapeFunction(Model,sym('X'));
Ha = Model.Ba.'*StrainTensor(Model)*Model.Ba;
Ma = Model.Ba.'*MassTensor(Model)*Model.Ba;
Model.Kee = double(int(Model.Phi.'*Ha*Model.Phi,sym('X'),0,1));
Model.Mee = double(int(Model.Phi.'*Ma*Model.Phi,sym('X'),0,1));
Model.Dee = Model.mu*Model.Kee;
Model.Phi = matlabFunction(Model.Phi);
Model.Nq = Model.NDof*Model.NModal;
end

%---------------------------------------------------------------------- set
function P = ShapeFunction(Model,X)

for ii = 1:Model.NModal, P(1,ii) = X.^(ii-1); end
for ii = 1:Model.NDof, Pc{ii,1} = P; end

P = blkdiag(Pc{:}) + 1e-16*X;
end

%---------------------------------------------------------------------- set
function dx = KinematicODE(Model,t,x)
Model.q = x;
Model.dq = 0*x;
Model.ddq = 0*x;

Model.t = t;
[~,Qc] = InverseDynamicModel(Model);
dx(:,1) = Model.Dee\(-Qc - Model.Kee*x);
end

%---------------------------------------------------------------------- set
function dx = DynamicODE(Model,t,x)
N = Model.NDof;
Model.q(:,1) = x(1:Model.NModal*Model.NDof);
Model.dq(:,1) = x(Model.NModal*Model.NDof+1:end);
Model.ddq(:,1) = Model.q(:,1)*0;
Model.t = t;

[~,Qc] = InverseDynamicModel(Model);

dx = x*0;
dx(1:Model.NModal*Model.NDof,1) = x(Model.NModal*Model.NDof+1:end);
dx(Model.NModal*Model.NDof+1:end,1) = Model.Mee\(-Model.Dee*Model.dq(:,1) ...
    -Qc - Model.Kee*Model.q(:,1));

end

%---------------------------------------------------------------------- set
function [F0,Q0] = InverseDynamicModel(Model)
%X = 0:0.01:1;
options = odeset('AbsTol',1e-3);

X = [0 1];

% forward integration
g0 = [1,0,0,0,0,0,0];
eta0 = [0,0,0,0,0,0];
deta0 = [0,0,0,0,0,0];

[~,yf] = ode45(@(t,x) ForwardODE(Model,t,x),X,...
    [g0 eta0 deta0],options);

% backwards integration
g1 = yf(end,1:7);
eta1 = yf(end,8:13);
deta1 = yf(end,14:19);
F1 = [0,0,1.0e-5,0,0,0];
Qa = zeros(Model.Nq,1).';

[~,yb] = ode45(@(t,x) BackwardODE(Model,t,x),flip(X),...
    [g1 eta1 deta1 F1 Qa],options);

F0 = -yb(end,20:25).';
Q0 = yb(end,26:end).';
end

%---------------------------------------------------------------------- set
function dg = ForwardODE(Model,t,g)

ee = Model.Ba*Model.Phi(t)*Model.q + Model.xia0;
dee = Model.Ba*Model.Phi(t)*Model.dq;
ddee = Model.Ba*Model.Phi(t)*Model.ddq;

Kappa = ee(1:3);
Gamma = ee(4:6);

Q = g(1:4);
r = g(5:7);
eta = g(8:13);
deta = g(14:19);

R = Quat2Rot(Q);
A = StrainMap(R*Kappa(:));

dg = zeros(19,1);
dg(1:4) = ((2*norm(Q))^(-1))*A*Q;
dg(5:7) = R*Gamma(:);
dg(8:13) = -admap(eta)*eta+dee;
dg(14:19) = -admap(eta)*deta - admap(deta)*eta + ddee;

end

%---------------------------------------------------------------------- set
function dg = BackwardODE(Model,t,g)

FBar = zeros(6,1);
ee = Model.Ba*Model.Phi(t)*Model.q + Model.xia0;
dee = Model.Ba*Model.Phi(t)*Model.dq;
ddee = Model.Ba*Model.Phi(t)*Model.ddq;

Q = g(1:4);
r = g(5:7);
eta = g(8:13)*0;
deta = g(14:19)*0;
L = g(20:25);
Kappa = ee(1:3);
Gamma = ee(4:6);

R = Quat2Rot(Q);
A = StrainMap(R*Kappa(:));
M = MassTensor(Model);
% 
% win = 0.1;
% 
% X0 = tdelta(t+win,win); 
% X0_ = tdelta(t-win,win);
% X1 = tdelta(t-1/2+win,win); 
% X1_ = tdelta(t-1/2-win,win);
% X2 = tdelta(t-1+win,win); 
% X2_ = tdelta(t-1-win,win);

% m1x = Model.P1(1)*(X1 - X0_);%*(smoothstep(f*Model.t) - smoothstep(f*Model.t-3));
% m2x = Model.P2(1)*(X2 - X1_);%*smoothstep(f*Model.t);
% m1y = Model.P1(2)*(X1 - X0_)*(smoothstep(f*Model.t) - ...
%     3*smoothstep(f*(Model.t-4.5)) + 2*smoothstep(f*(Model.t-7.5)));
% m2y = Model.P2(2)*(X2 - X1_)*(smoothstep(f*Model.t) - ...
%     2*smoothstep(f*(Model.t-4.2)) + smoothstep(f*(Model.t-7.2)));

% m1x = Model.P1(1)*(X1 - X0_);
% m2x = Model.P2(1)*(X2 - X1_);
% m1y = Model.P1(2)*(X1 - X0_);
% m2y = Model.P2(2)*(X2 - X1_);
% 
p0 = Model.Density;
A = pi*Model.Radius^2;
% 
%FBar = 0*[0,m1x+m2x,m1y+m2y,0,0,0].';

%FBar(4:end) = -r(1)*p0*A*Model.Grav*[1;0;0];
%Fc = 0*[cross(D1,Gamma_c); Gamma_c]*(a*ndelta(t-0.5,0.01));

dg = 0*g;
dg(1:4) = ((2*norm(Q))^(-1))*A*Q;
dg(5:7) = R*Gamma(:);
dg(8:13) = -admap(eta)*eta+dee;
dg(14:19) = -admap(eta)*deta - admap(deta)*eta + ddee;
dg(20:25) = (admap(ee)).'*L + M*deta - admap(eta).'*M*eta - FBar;
dg(26:end) = -Model.Phi(t).'*Model.Ba.'*L;

end
    
end
end

%---------------------------------------------------------------------- set
function C = HatOperator(a)
C = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
end

%---------------------------------------------------------------------- set
function g = admap(x)
W = x(1:3); U = x(4:6); g = zeros(6);
Wh = HatOperator(W); Uh = HatOperator(U);
g(1:3,1:3) = Wh;
g(4:6,4:6) = Wh;
g(4:6,1:3) = Uh;
end

%---------------------------------------------------------------------- set
function Adg = Adgmap(Q,r)
R = Quat2Rot(Q);
Adg = zeros(6); 
Adg(1:3,1:3) = R;
Adg(4:6,4:6) = R;
Adg(4:6,1:3) = HatOperator(r)*R;
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

function H = StrainTensor(Model)
A = pi*Model.Radius^2;
J1 = 0.5*pi*Model.Radius^4;
J2 = 0.25*pi*Model.Radius^4;
J3 = 0.25*pi*Model.Radius^4;
E0 = Model.E;
nu0 = Model.nu;
G0 = (E0)/(2*(1+nu0));
H = diag([G0*J1,E0*J2,E0*J3,E0*A, G0*A, G0*A]);
end

function M = MassTensor(Model)
p0 = Model.Density;
A = pi*Model.Radius^2;
J1 = 0.5*pi*Model.Radius^4;
J2 = 0.25*pi*Model.Radius^4;
J3 = 0.25*pi*Model.Radius^4;
M = diag([p0*J1,p0*J2,p0*J3,p0*A, p0*A, p0*A]);
end