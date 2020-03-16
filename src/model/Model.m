classdef Model

    properties (Access = public)
        Links;
        Dof;
        NDof;
        NLink = 1;
        tspan = 1;
        Table;
        g;
        Controller;
    end
    
    properties (Access = private)
        Nq;
        q; 
        dq;
        ddq;
        t;
        E;
        nu = 0.45;
        mu = 0.75;
        xia0 = [0,0,0,1,0,0].';
        Phi;
        Ba; Bc; 
        Mee; Kee; Dee; Qee;
        NModal;
        Grav;
        Radius;
        Density;
        Length ;
        Length0;
        P1; P2; P3;
        LumpedMass;
        Chebyshev;
        DistLoad;
        PressureArea;
        
        Movie;
        MovieStart;
        MovieAxis;
    end
    
%--------------------------------------------------------------------------
methods  
    
%-------------------------------------------------------------- Model Class
function obj = Model(Table,varargin) 
    
    obj.Table = Table;
    obj.NLink = size(Table,1);
    obj.LumpedMass = false;
    obj.Grav = 9.81;
    obj.Movie = false;
    obj.MovieStart = false;
    obj.MovieAxis = 1.5*[-.5 .5 -.5 .5 -1 0.5];
    obj.Density = 0.01;
    
    obj.Radius = 1e-2;
    obj.E = 4e2;
    obj.nu = 0.4;
    obj.mu = 0.3;
    obj.Length = 1;
    obj.NModal = 2;
    obj.Chebyshev = true;
    obj.PressureArea = 5e-5;
    
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
function Model = generate(Model)
    Model.Length0 = Model.Length;
    Model = GenerateCosserat(Model);
end

%---------------------------------------------------------------------- set
function Model = solve(Model)
    
dt = 1e-1;
T = 0:dt:Model.tspan;
q0 = zeros(Model.NModal*Model.NDof,1);

opts = optimoptions('fsolve','Display','none','TolFun',1e-5);
ys = fsolve(@(x) KinematicODE(Model,T,x),q0,opts);

Model.g = ys(:)';
Model.t = 1;

end

%----------------------------------------------------------------- simulate
function Model = simulate(Model)
N = Model.NModal*Model.NDof;
dt = 1e-1;
T = 0:dt:Model.tspan;
q0 = zeros(N,1);
dq0 = zeros(N,1);

opt = odeset('RelTol',1e-1,'AbsTol',1e-1);
[ts,ys] = ode23tb(@(t,x) DynamicODE(Model,t,x),T,[q0 dq0],opt);

Model.g = ys;
Model.t = ts;
end

%--------------------------------------------------------------------- show
function show(Model)
   
N = Model.Nq;
X = linspace(0,1,200);
FPS = round((1/12)/(mean(diff(Model.t))));
axis([-.5 .5 -.5 .5 -1 .5]);

if length(Model.t) > 2
FPS = round((1/12)/(mean(diff(Model.t))));
else, FPS = 1;
end

for ii = 1:FPS:length(Model.t)
    Model.q(:,1) = Model.g(ii,1:N).';
    Model.dq(:,1) = 0*Model.q(:,1);
    Model.ddq(:,1) = 0*Model.q(:,1);
    ee = Model.Ba*Model.Phi(Model.Length)*Model.q;
    Model.Length = ee(4);
    g0 = [1,0,0,0,0,0,0];
    eta0 = [0,0,0,0,0,0];
    deta0 = [0,0,0,0,0,0];
    [~,yf] = ode45(@(t,x) ForwardODE(Model,t,x,Model.q(:,1),Model.q(:,1)*0,...
        Model.q(:,1)*0),X,[g0 eta0 deta0]);
    G = yf(:,1:7);
    
    figure(101);
    cla;
    plot3(G(:,6),G(:,7),-G(:,5),'linewidth',2,'Color',col(1));
    axis equal
    
end
    
end

%------------------------------------------------------------ show 3D model
function showModel(Model)
    
figure(101); hold all; cla;
N = Model.Nq;
Model.Length = Model.Length0;
X = linspace(0,1,200);

%msh = Gmodel('SlenderRod.stl');
msh = Gmodel('SoftActuatorRedux.stl');
mshgr = Gmodel('SoftGripperRedux.stl');
%mshgr = Gmodel([]);
%mshgr = mshgr.set('Node0',mshgr.Node*1e-6);
% mshgr = mshgr.set('Node',mshgr.Node*1e-6);

% set texture
msh.Texture = 1.25*grey;%1.25*base;
mshgr.Texture = 1.25*grey;%1.25*base;

msh = msh.bake();
mshgr = mshgr.bake();
msh = msh.render();
mshgr = mshgr.render();

LinkID = knnsearch(X(:),msh.Node(:,3));
axis equal;
if ~isempty(Model.MovieAxis), axis(Model.MovieAxis); end
drawnow;
%axis([-.5 .5 -.5 .5 -1 .5])
view(60,5);
msh.update();
mshgr.update();

if length(Model.t) > 2, FPS = round((1/12)/(mean(diff(Model.t)))); i0 = 1;
else, FPS = 1; i0 = 1;
end

for ii = i0:FPS:length(Model.t)

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
    
    [~,yf] = ode45(@(t,x) ForwardODE(Model,t,x,Model.q(:,1),Model.q(:,1)*0,...
        Model.q(:,1)*0),X,[g0 eta0 deta0]);

    SweepSE3 = yf(:,1:7);
    
    msh = Blender(msh,'Rotate',{'z',30});
    mshgr = Blender(mshgr,'Rotate',{'z',30});
    msh = Blender(msh,'Sweep', {LinkID,SweepSE3});
    msh = Blender(msh,'Scale',{'z',-1});
    mshgr = Blender(mshgr,'SE3',yf(end,1:7));
    mshgr = Blender(mshgr,'Scale',{'z',-1});
    
    mshgr.updateNode();
    msh.updateNode();
    msh.update();
    mshgr.update();
        
    if ~isempty(Model.MovieAxis), axis(Model.MovieAxis); end
    
    if Model.Movie
        if Model.MovieStart == false
            Model.MovieStart = true;
            MovieMaker(Model,'mdl','Start');
        else
            MovieMaker(Model);
        end
    end
    
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
N = Model.NModal*Model.NDof;
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

if ~Model.Chebyshev
for ii = 1:Model.NModal, P(1,ii) = X.^(ii-1); end
else
for ii = 1:Model.NModal, P(1,ii) = chebyshev(X,ii-1); end
end
for ii = 1:Model.NDof, Pc{ii,1} = P; end

P = blkdiag(Pc{:}) + 1e-16*X;
end

%---------------------------------------------------------------------- set
function P = RectDelta(~,x,h)
P = piecewise(abs(x) < h/2, 1/h, 1);
end

%---------------------------------------------------------------------- set
function P = Dirac(~,x,n)
P = dirac(x,n);
end

%---------------------------------------------------------------------- set
function dx = KinematicODE(Model,~,x)
N = Model.NModal*Model.NDof;
Model.t = Model.tspan;
Model.q = x;
q_   = x(1:N); 
g0   = [1,0,0,0,0,0,0];
n0   = [0,0,0,0,0,0];

% static potential/reaction forces
[~,Qc,yf,ts] = InverseDynamicModel(Model,q_,0*q_,0*q_,g0,n0,n0);

[J0,J1] = JacobianCompute(Model,ts, yf(:,1:7));

t1 = PressureVector(1,0,1);
t2 = PressureVector(0,0,1);
Qs = (J1.'*(t1 - t2) + J0.'*t2);

Qa = Qc + Qs;

dx(:,1) = Model.Dee\(-Qa - Model.Kee*q_);
end

%---------------------------------------------------------------------- set
function dx = DynamicODE(Model,t,x)
N = Model.NModal*Model.NDof;
Model.q(:,1)  = x(1:N);
Model.dq(:,1) = x(N+1:end);
Model.t = t;

q_   = x(1:N);
dq_  = x(N+1:end);
ddq_ = Model.dq(:,1)*0;
g0   = [1,0,0,0,0,0,0];
n0   = [0,0,0,0,0,0];

% static potential/reaction forces
[~,Qc] = InverseDynamicModel(Model,q_,0*dq_,0*ddq_,g0,n0,n0);

% dynamic coriolis/viscous forces
[~,Qr,yf,ts] = InverseDynamicModel(Model,q_,dq_,0*ddq_,g0,n0,n0);
Qv = Qr - Qc;

% building the interia matrix
MC = zeros(N);
for ii = 1:(Model.NModal*Model.NDof)
    [~,Mtmp] = InverseDynamicModel(Model,q_,0*dq_,dalpha(ii,N),g0,n0,n0);
    MC(:,ii) = -Mtmp + Qc;
end

% construct jacobian matrix
[J0,J1] = JacobianCompute(Model,ts, yf(:,1:7));

% compute control action
g1 = yf(end,1:7);
eta1 = yf(end,8:13);
Qd = -J0.'*Model.Controller(Model.t,g1);

% compute pressure inputs
H = PressureMatrix;
I = ones(6,1);
O = zeros(6,1);
opt = optimoptions('fmincon','Algorithm','sqp','Display','none');
tau = fmincon(@(x) x.'*x,I*1e-3,[J1.'*H, (J0.'- J1.')*H],Qd,[],[],O,I,[],opt);
if isempty(tau), tau = zeros(6,1); end

t1 = H*tau(1:3);
t2 = H*tau(4:6);
Qs = (J1.'*(t1 - t2) + J0.'*t2);

Qa = Qc + Qv + Qs;

dx = x;
dx(1:N,1) = x(Model.NModal*Model.NDof+1:end);
dx(N+1:end,1) = MC\(-Qa - Model.Kee*q_ -Model.Dee*dq_ );

% delta operator
function y = dalpha(a,N), y = logical(1:N==a).'; end
    
end

%---------------------------------------------------------------------- set
function [F0,Q0,Y,X] = InverseDynamicModel(Model,q,dq,ddq,g0,eta0,deta0)
X = linspace(0,Model.Length,30);

% forward integration
[~,yf] = oderk(@(t,x) ForwardODE(Model,t,x,q,dq,ddq),X,[g0, eta0 deta0]);

g1    = yf(end,1:7);
eta1  = yf(end,8:13);
deta1 = yf(end,14:19);

F1 = zeros(1,6);
Qa = zeros(Model.Nq,1).';

% backwards integration
[~,yb] = oderk(@(t,x) BackwardODE(Model,t,x,q,dq,ddq),flip(X),...
     [g1,eta1,deta1,F1,Qa]);

Y = flipud(yb(:,1:19));
F0 = -yb(end,20:25).';
Q0 = yb(end,26:end).';
end

%--------------------------------------- forwards integration of kinematics
function dg = ForwardODE(Model,t,g,q,dq,ddq)
ee = Model.Ba*Model.Phi(t)*q + Model.xia0;
dee = Model.Ba*Model.Phi(t)*dq;
ddee = Model.Ba*Model.Phi(t)*ddq;

Kappa = ee(1:3);
Gamma = ee(4:6);
Q     = g(1:4);
r     = g(5:7);
eta   = g(8:13);
deta  = g(14:19);

R = Quat2Rot(Q);
A = StrainMap(R*Kappa(:));

dg = zeros(19,1);

dg(1:4)   = ((2*norm(Q))^(-1))*A*Q;
dg(5:7)   = R*Gamma(:);
dg(8:13)  = -admap(ee)*eta+dee;
dg(14:19) = -admap(ee)*deta - admap(dee)*eta + ddee;

end

%----------------------------------------- backward integration of dynamics
function dg = BackwardODE(Model,t,g,q,dq,ddq)
ee   = Model.Ba*Model.Phi(t)*q + Model.xia0;
dee  = Model.Ba*Model.Phi(t)*dq;
ddee = Model.Ba*Model.Phi(t)*ddq;

Q      = g(1:4);
r      = g(5:7);
eta    = g(8:13);
deta   = g(14:19);
Lambda = g(20:25);
Kappa  = ee(1:3);
Gamma  = ee(4:6);

R    = Quat2Rot(Q);
A    = StrainMap(R*Kappa(:));
M    = MassTensor(Model);
adee = admap(ee);
adet = admap(eta);

dg = 0*g;
FBar = zeros(6,1);
FBar(4:6) = FBar(4:6) - R.'*M(4,4)*[Model.Grav,0,0].';

dg(1:4)    = ((2*norm(Q))^(-1))*A*Q;
dg(5:7)    = R*Gamma(:);
dg(8:13)   = -adee*eta + dee;
dg(14:19)  = -adee*deta - admap(dee)*eta + ddee;
dg(20:25)  = (adee).'*Lambda + M*deta - adet.'*M*eta - FBar;
dg(26:end) = -Model.Phi(t).'*Model.Ba.'*Lambda;

end

%-------------------------------------------------- compute jacobian matrix
function [J0, J1] = JacobianCompute(Model,t,gmat)
dx = mean(diff(t));

G0 = gmat(1,:);
GL = gmat(end,:);

J0 = dx*0.5*localJ(Model,G0,0);

for ii = 1:length(t)-1
    J0 = J0 + dx*localJ(Model,gmat(ii,:),t(ii));
    if ii == floor(length(t)/2)-1
    GMID = gmat(ii+1,:);
    J1 = J0 + 0.5*dx*localJ(Model,GMID,t(ii+1));
    J1 = transpose(Ad(GMID))*J1;
    end
end

J0 = transpose(Ad(GL))*(J0 + 0.5*dx*localJ(Model,GL,1));

%----------------------------------------------------------------
function J = localJ(Model,g,t), Q = g(1:4); r = g(5:7);
        J = Adgmap(Q,r)*Model.Ba*Model.Phi(t);
end
%----------------------------------------------------------------
function A = Ad(g), Q = g(1:4); r = g(5:7); A = Adgmap(Q,r); end
%----------------------------------------------------------------
    
end
  
end
end

%----------------------------------------- ODE solver // runge kutta method
function [t,dx] = oderk(f,t,y0)
N = length(t);
F = @(x,y) f(x,y);
y = zeros(length(y0),N);
y(:,1) = y0(:);
h = mean(diff(t));

for i = 1:(N-1)  
    ti = t(i);
    yi = y(:,i);
    
    k1 = F(ti,yi);
    k2 = F(ti+0.5*h,yi+0.5*h*k1);
    k3 = F(ti+0.5*h,yi+0.5*h*k2);
    k4 = F(ti+h,yi+k3*h);
    
    y(:,i+1) = yi + (1/6)*(k1+2*k2+2*k3+k4)*h;  
end

dx = transpose(y);
t = t(:);
end

%--------------------------------------- ODE solver // forward euler method
function [t,dx] = ode1(f,t,y0)
N = length(t);
F = @(x,y) f(x,y);
y = zeros(length(y0),N);
y(:,1) = y0(:);
h = mean(diff(t));

for i = 1:(N-1)  
    k1 = F(t(i),y(:,i));
    y(:,i+1) = y(:,i) + k1*h;  
end

dx = transpose(y);
t = t(:);
end

%---------------------------------------- isomorhphism between SO(3) and R6
function y = isomSO3(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; x3, 0, -x1; -x2, x1, 0];
end

%----------------------------------------- adjoint action on lie alg. se(3)
function g = admap(x)
W = [x(1);x(2);x(3)];
U = x(4:6); 
g = zeros(6);
Wh = isomSO3(W); 
Uh = isomSO3(U);
g(1:3,1:3) = Wh;
g(4:6,4:6) = Wh;
g(4:6,1:3) = Uh;
end

%------------------------------------------- adjoint map on lie group SE(3)
function Adg = Adgmap(Q,r)
R = Quat2Rot(Q);
Adg = zeros(6); 
Adg(1:3,1:3) = R;
Adg(4:6,4:6) = R;
Adg(4:6,1:3) = isomSO3(r)*R;
end

%----------------------------------------------------------- strain mapping
function A = StrainMap(K)
k1 = K(1); k2 = K(2); k3 = K(3);
A = [ 0, -k1, -k2, -k3; k1,   0, -k3,  k2; 
     k2,  k3,   0, -k1; k3, -k2,  k1,  0];
end

%----------------------------------------------- quaterion to rotation mat.
function R = Quat2Rot(q)
w = q(1); x = q(2); y = q(3); z = q(4);
Rxx = 1 - 2*(y^2 + z^2); Rxy = 2*(x*y - z*w); Rxz = 2*(x*z + y*w); 
Ryx = 2*(x*y + z*w); Ryy = 1 - 2*(x^2 + z^2); Ryz = 2*(y*z - x*w );
Rzx = 2*(x*z - y*w ); Rzy = 2*(y*z + x*w ); Rzz = 1 - 2 *(x^2 + y^2);

R = [Rxx, Rxy, Rxz; Ryx, Ryy, Ryz; Rzx, Rzy, Rzz];
end

%----------------------------------------- infinitesimal body strain tensor 
function H = StrainTensor(Model)
A = pi*Model.Radius^2;
J1 = 0.5*pi*Model.Radius^4;
J2 = 0.25*pi*Model.Radius^4;
J3 = 0.25*pi*Model.Radius^4;
E0 = Model.E;
nu0 = Model.nu;
G0 = (E0)/(2*(1+nu0));
H = diag([G0*J1,E0*J2,E0*J3, 0.001*E0*A, G0*A, G0*A]);
end

%---------------------------------------- infinitesimal body inertia tensor 
function M = MassTensor(Model)
p0 = Model.Density;
A = pi*Model.Radius^2;
J1 = 0.5*pi*Model.Radius^4;
J2 = 0.25*pi*Model.Radius^4;
J3 = 0.25*pi*Model.Radius^4;
M = diag([p0*J1,p0*J2,p0*J3,p0*A, p0*A, p0*A]);
end

%----------------------------------------------------- control input tensor 
function P = PressureMatrix
%P1 = [p11;p12;p13];
A = 5e-6;
r = 1;

H = [ 0,              0,             0;
      0, -0.5*r*sqrt(3), 0.5*r*sqrt(3);
     -r,          0.5*r,         0.5*r;
      0,              0,             0;
      0,              0,             0;
      0,              0,             0];

P = A*H;
  
end

%-------------------------------------------------------------- movie maker
function MovieMaker(Model,Name,Request)
if nargin < 2, Request = ''; end

if Model.Movie
    switch(Request)
        case('Start')
            filename = string([Name,'_', char(datetime(now,...
                'ConvertFrom','datenum')),'.gif']);
            
            filename = erase(filename,[":"," "]);
            background(metropolis);
            if ~isempty(Model.MovieAxis), axis(Model.MovieAxis); end
            drawnow;
            gif(char(filename),'frame',gcf,'nodither');
        otherwise
            background(metropolis);
            if ~isempty(Model.MovieAxis), axis(Model.MovieAxis); end
            drawnow;
            gif;
    end
end
end
