classdef Model

    properties (Access = public)
        Dof;
        NDof;
        NModal;
        NDisc;
        Table;
        g;
        ge;
        tau;
        Texture;
        q;
        dq;
        q0;
        dq0;
        Gain;
        Point;
        Controller;
        Area, Jxx, Jyy, Jzz;
    end
    
    properties (Access = private)
        Nq;
        xd;
        
        t;
        E, Nu, Mu;
        xia0 = [0,0,0,1,0,0].';
        Phi;
        Ba; Bc; 
        SpaceStep;
        TimeStep
        Gravity;
        Radius;
        Density;
        Tdomain;
        Sdomain
        Length;
        Length0;
        SimTime;
        DistLoad;
        PArea;
        
        Movie;
        MovieStart;
        MovieAxis;
    end
    
%--------------------------------------------------------------------------
methods  
    
%-------------------------------------------------------------- Model Class
function obj = Model(Table,varargin) 
    
    obj.Table  = Table;
    obj.Controller = 1;
    %obj.NLink = size(Table,1);
    obj.Movie      = false;
    obj.MovieStart = false;
    obj.MovieAxis  = [];
    
    obj.NModal    = 2;
    obj.NDisc     = 2;
    obj.SpaceStep = 11;
    obj.TimeStep  = 0.01;
    obj.Tdomain   = 25;
    obj.Sdomain   = 1;
    
    obj.Density = 0.01;
    obj.Radius  = 0.01;
    obj.E       = 1e5;
    obj.Nu      = 0.4;
    obj.Mu      = 0.1;
    obj.Length  = 1;
    obj.Gravity = 0;
    obj.PArea   = 5e-5;
    obj.Area    = pi*obj.Radius^2;
    obj.Jxx     = 0.5*obj.Density*obj.Radius^2;
    obj.Jyy     = 0.25*obj.Density*obj.Radius^2;
    obj.Jzz     = 0.25*obj.Density*obj.Radius^2;
    
    obj.Point = [1,0,0,0,1,0,0];
    obj.Gain = [1e-7,0.15];

    obj.Texture = prusa;
    
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
    Model = GenerateConfigFile(Model);
end

%---------------------------------------------------------------------- set
function Model = solve(Model)   
x0 = zeros(Model.NModal*Model.NDof,1);
dt = 1e-1;
opts = optimoptions('fsolve','Display','none','TolFun',1e-3,...
    'Algorithm','levenberg-marquardt');
ys = fsolve(@(x) KinematicODE(Model,Model.tspan,x),x0,opts);

% opt = odeset('RelTol',1e-1,'AbsTol',1e-1);
% [~,ys] = oderk(@(t,x) KinematicODE(Model,Model.tspan,x),T,q0);

Model.g = ys(:).';
Model.t = 1;
end

%--------------------------------------------------------------------- make
function Model = make(Model)
tic
CMD = 'cd src/model/tools/solver/build && make';
system(CMD);
toc    
end

%---------------------------------------------------------------------- set
function Model = csolve(Model)   
dir_path = './src/model/tools/solver/build';

%//////////////////////////////////
out = fullfile(dir_path,'log/state.log');
fileID = fopen(out,'w');
fprintf(fileID,'%d\n',Model.q0);
fclose(fileID);
%//////////////////////////////////
out = fullfile(dir_path,'log/state_dt.log');
fileID = fopen(out,'w');
fprintf(fileID,'%d\n',Model.dq0);
fclose(fileID);

tic
CMD = 'cd src/model/tools/solver/build && ./solver config.m';
system(CMD);
toc
y = load('src/model/tools/solver/build/log/state.log');
Model.q = y(:,2:end);
Model.t = y(:,1);
y = load('src/model/tools/solver/build/log/tau.log');
Model.tau = y(:,2:end);
y = load('src/model/tools/solver/build/log/g.log');
Model.ge = y(:,2:end);
y = load('src/model/tools/solver/build/log/statelog_dt.log');
Model.dq = y(:,2:end);
y = load('src/model/tools/solver/build/log/xd.log');
Model.xd = y(:,2:end);

end

%----------------------------------------------------------------- simulate
function Model = simulate(Model)
N = Model.NModal*Model.NDof;
dt = 1e-2;
T = 0:dt:Model.tspan;
x0 = zeros(N,1);
dx0 = zeros(N,1);

x0(1) = 0;

%if ~isempty(Model.q0), x0 = Model.q0; end

opt = odeset('RelTol',1e-1,'AbsTol',1e-1);
[ts,ys] = ode23tb(@(t,x) DynamicODE(Model,t,x),T,[x0;dx0],opt);
%[ts,ys] = ode1(@(t,x) DynamicODE(Model,t,x),T,[q0; dq0]);

Model.g = ys;
Model.t = ts;
end

%----------------------------------------------------------------- simulate
function Model = inertia(Model)
    Model.Area    = pi*Model.Radius^2;
    Model.Jxx     = 0.5*Model.Density*Model.Radius^2;
    Model.Jyy     = 0.25*Model.Density*Model.Radius^2;
    Model.Jzz     = 0.25*Model.Density*Model.Radius^2; 
end

%--------------------------------------------------------------------- show
function show(Model)
N = Model.Nq;
X = linspace(0,1,20);

if length(Model.t) > 2, FPS = round((1/13)/(mean(diff(Model.t))));
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
    if ~isempty(Model.MovieAxis), axis equal; axis(Model.MovieAxis); end
end
    
end

%------------------------------------------------------------ show 3D model
function showModel(Model)
    
figure(101); 
hold all; 
N = Model.Nq;
%Model.Length = Model.Length0;
X = linspace(0,Model.Sdomain,200);

%msh = Gmodel('SlenderRod.stl');
%msh = Gmodel('SoftActuatorRedux.stl');
%msh = Gmodel('SoftActuatorPlanarRedux.stl');
msh = Gmodel('Pneulink.stl'); 
%msh = Gmodel('Pneunet.stl'); 
assignin('base','msh',msh);
mshgr = Gmodel('SoftGripperRedux.stl'); assignin('base','mshgr',mshgr);
%mshgr = Gmodel([]);
%msh = msh.set('Node0',mshgr.Node);
mshgr = mshgr.set('Node0',mshgr.Node*1.0e-5);
mshgr = mshgr.set('Node',mshgr.Node*1.0e-5);
msh = msh.set('Node0',msh.Node*0.064);
msh = msh.set('Node',msh.Node*0.064);

% set texture
msh.Texture = Model.Texture;
mshgr.Texture = Model.Texture;

msh = msh.bake();
mshgr = mshgr.bake();
msh   = msh.render();
mshgr = mshgr.render();
axis equal; 

LinkID = knnsearch(X(:),msh.Node(:,3));
if isempty(Model.MovieAxis)
    Model.MovieAxis = boxhull(msh.Node); 
end
axis(Model.MovieAxis);
drawnow;
%axis([-.5 .5 -.5 .5 -1 .5])
view(50,10);
msh.update();
mshgr.update();
msh.ground(Model.MovieAxis);

%background('k');

if length(Model.t) > 2, FPS = max(round((1/12)/(mean(diff(Model.t)))),1); i0 = 1;
else, FPS = 1; i0 = 1;
end

plotvector([0;0;0],[0.2;0;0],'Color',col(1),'linewidth',2,'MaxHeadSize',0.75);
plotvector([0;0;0],[0;0.2;0],'Color',col(2),'linewidth',2,'MaxHeadSize',0.75);
plotvector([0;0;0],[0;0;0.2],'Color',col(3),'linewidth',2,'MaxHeadSize',0.75);

if ~isempty(Model.xd)
% p = Model.point(4:6).';
p = Model.xd(end,5:7);
R = quat2rot(Model.xd(end,1:4));
fr = [p(1,3);p(1,2);-p(1,1)];
% % 
f0 = plotvector(fr,R*[0.2;0;0],'Color',col(1),'linewidth',2,'MaxHeadSize',0.75);
f1 = plotvector(fr,R*[0;0.2;0],'Color',col(2),'linewidth',2,'MaxHeadSize',0.75);
f2 = plotvector(fr,R*[0;0;0.2],'Color',col(3),'linewidth',2,'MaxHeadSize',0.75);
end

vline = [];

for ii = [1:FPS:length(Model.t),length(Model.t)]

    msh.reset();
    mshgr.reset();

    Qtmp = Model.q(ii,1:N).';
    ee = Model.Ba*Model.Phi(Model.Length)*Qtmp;
    
    Model.Length = ee(4);
    
    g0 = [1,0,0,0,0,0,0];
    eta0 = [0,0,0,0,0,0];
    deta0 = [0,0,0,0,0,0];
    
    [~,yf] = ode23(@(t,x) ForwardODE(Model,t,x,Qtmp,...
         Qtmp*0, Qtmp*0),X,[g0 eta0 deta0]);

    SweepSE3 = yf(:,1:7);
    
    msh = Blender(msh,'Rotate',{'z',-30});
    %msh = Blender(msh,'Scale',0.064);
    mshgr = Blender(mshgr,'Rotate',{'z',-30});
    msh = Blender(msh,'Sweep', {LinkID,SweepSE3});
    msh = Blender(msh,'Scale',{'z',-1});
    %msh = Blender(msh,'Rotate',{'y',90});
    %msh = Blender(msh,'Rotate',{'z',180});
    mshgr = Blender(mshgr,'SE3',yf(end,1:7));
    mshgr = Blender(mshgr,'Scale',{'z',-1});
    
    %vline = [vline; SweepSE3(end,7),SweepSE3(end,6),SweepSE3(end,5),Model.t(ii)];

    if ii == 1      
       h = plot3(SweepSE3(:,7),SweepSE3(:,6),-SweepSE3(:,5),'linewidth',1.5,...
       'Color',col(1));
       hm = plot3(SweepSE3(end,7),SweepSE3(end,6),-SweepSE3(end,5),'.',...
           'markersize',10,'Color',col(1));
       
    if ~isempty(Model.xd)
       p = Model.xd(:,4:6);
       h0 = plot3(p(:,3),p(:,2),-p(:,1),'k--');
    end
        
    else
        
        delete(h);
         if ~isempty(Model.xd)
        delete(hm);
        delete(h0);
        end
       if size(vline,1) > 4, delete(hvm); end
       
        
        if ~isempty(Model.xd)
        delete(h0);
        delete(f0);
        delete(f1);
        delete(f2);
        %delete(h1);
        %delete(h2);
        p = Model.xd(ii,5:7);
        R = quat2rot(Model.xd(ii,1:4));
        fr = [p(3);p(2);-p(1)];

         f0 = plotvector(fr,R*[0.2;0;0],'Color',col(1),'linewidth',2,'MaxHeadSize',0.75);
         f1 = plotvector(fr,R*[0;0.2;0],'Color',col(2),'linewidth',2,'MaxHeadSize',0.75);
         f2 = plotvector(fr,R*[0;0;0.2],'Color',col(3),'linewidth',2,'MaxHeadSize',0.75);
        end
        
        h = plot3(SweepSE3(:,7),SweepSE3(:,6),-SweepSE3(:,5),'linewidth',2,...
        'Color',col(2));
    
        %hm = plot3(SweepSE3(end,7),SweepSE3(end,6),-SweepSE3(end,5),'.',...
        %   'markersize',10,'Color',col(1));
       
        %h0 = plot3(p(:,3),p(:,2),-p(:,1),'k--');

    end
      
    %mshgr.updateNode();
    %msh.updateNode();
    msh.update();
    mshgr.update();
    
    title(['T = ',num2str(Model.t(ii),3)]);
    %title(['iterations = ',num2str(ii,3)]);
        
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

Model.Ba   = I6(:,xa);
Model.Bc   = I6(:,xc);
Model.NDof = size(Model.Ba,2);
Model.Phi  = @(x) ShapeFunction(Model,x);
Model.Nq   = Model.NDof*Model.NModal;
Model.q0   = zeros(Model.NDof*Model.NModal,1);
Model.dq0  = zeros(Model.NDof*Model.NModal,1);
end

%---------------------------------------------------------------------- set
function Model = GenerateConfigFile(Model)
   
global FID;
File = [cdsoro,'/src/model/tools/solver/build/config.m'];
delete(File);
FID = fopen(File,'w');

fprintf(FID,'[options] \n');
if(Model.Controller)
    fprintf(FID,'KINEMATIC_CONTROLLER = 0 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 1 \n');
else
    fprintf(FID,'KINEMATIC_CONTROLLER = 0 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 0 \n');
end
fprintf(FID,'WRITE_OUTPUT         = 1 \n');
%fprintf(FID,'POINT_INPUT          = 0 \n');

if(Model.NDisc > 1)
    fprintf(FID,'DISCONTINIOUS        = 1 \n');
else
    fprintf(FID,'DISCONTINIOUS        = 0 \n');
end

fprintf(FID,'\n[cosserat] \n');
fprintf(FID,['K1 = ', num2str(Model.Table(1)), '\n']);
fprintf(FID,['K2 = ', num2str(Model.Table(2)), '\n']);
fprintf(FID,['K3 = ', num2str(Model.Table(3)), '\n']);
fprintf(FID,['E1 = ', num2str(Model.Table(4)), '\n']);
fprintf(FID,['E2 = ', num2str(Model.Table(5)), '\n']);
fprintf(FID,['E3 = ', num2str(Model.Table(6)), '\n']);

fprintf(FID,'\n[model] \n');
fprintf(FID,['NMODE   = ', num2str(Model.NModal), '\n']);
fprintf(FID,['NDISC   = ', num2str(Model.NDisc), '\n']);
fprintf(FID,['SDOMAIN = ', num2str(Model.Length), '\n']);
fprintf(FID,['TDOMAIN = ', num2str(Model.Tdomain), '\n']);

fprintf(FID,'\n[solver] \n');
fprintf(FID,['SPACESTEP = ', num2str(Model.SpaceStep), '\n']);
fprintf(FID,['TIMESTEP  = ', num2str(Model.TimeStep), '\n']);
fprintf(FID,['INTSTEP   = ', '100', '\n']);
fprintf(FID,['ATOL      = ', '1e-2', '\n']);
fprintf(FID,['RTOL      = ', '1e-2', '\n']);
fprintf(FID,['LAMBDA    = ', '1e-2', '\n']);
fprintf(FID,['MAX_IMPL  = ', '-1', '\n']);
fprintf(FID,['MAX_ITER  = ', '1', '\n']);
fprintf(FID,['MAX_IK    = ', '1500', '\n']);
fprintf(FID,['SPEEDUP   = ', '80', '\n']);

fprintf(FID,'\n[physics] \n');
fprintf(FID,['RHO      = ', num2str(Model.Density), '\n']);
fprintf(FID,['EMOD     = ', num2str(Model.E), '\n']);
fprintf(FID,['NU       = ', num2str(Model.Nu), '\n']);
fprintf(FID,['MU       = ', num2str(Model.Mu), '\n']);
fprintf(FID,['PRS_AREA = ', '1e-5', '\n']);
fprintf(FID,['GRAVITY  = ', num2str(Model.Gravity), '\n']);
fprintf(FID,['RADIUS   = ', '0.01', '\n']);
fprintf(FID,['AREA     = ', num2str(Model.Area), '\n']); %'3.0000e-04'
fprintf(FID,['J_XX     = ', num2str(Model.Jxx), '\n']); % '2.7250e-10'
fprintf(FID,['J_YY     = ', num2str(Model.Jyy), '\n']); % '2.2500e-11'
fprintf(FID,['J_ZZ     = ', num2str(Model.Jzz), '\n']); %'2.5000e-10'

fprintf(FID,'\n[control] \n');
fprintf(FID,['KP = ', num2str(Model.Gain(1)), '\n']);
fprintf(FID,['KD = ', num2str(Model.Gain(2)), '\n']);

fprintf(FID,'\n[setpoint] \n');
fprintf(FID,['Q1d = ', num2str(Model.Point(1)), '\n']);
fprintf(FID,['Q2d = ', num2str(Model.Point(2)), '\n']);
fprintf(FID,['Q3d = ', num2str(Model.Point(3)), '\n']);
fprintf(FID,['Q4d = ', num2str(Model.Point(4)), '\n']);
fprintf(FID,['Xd =  ', num2str(Model.Point(5)), '\n']);
fprintf(FID,['Yd =  ', num2str(Model.Point(6)), '\n']);
fprintf(FID,['Zd =  ', num2str(Model.Point(7)), '\n']);

fclose('all');
    
end

%---------------------------------------------------------------------- set
function P = ShapeFunction(Model,X)

Pc = cell(Model.NDof,1);
X = single(X);

if Model.NDisc>1
    P11 = zeros(1,Model.NModal/Model.NDisc);
    P22 = zeros(1,Model.NModal/Model.NDisc);
else
    P11 = zeros(1,Model.NModal);
end


if Model.NDisc>1
    for ii = 1:(Model.NModal/Model.NDisc)
         P11(1,ii) = sign(z1(X)).*legendre(z1(X),ii-1);
         P22(1,ii) = sign(z2(X)).*legendre(z2(X),ii-1);
         %P11(1,ii) = w1(X)*chebyshev(z1(X),ii-1);
         %P22(1,ii) = w2(X)*chebyshev(z2(X),ii-1);

    end
    P = horzcat(P11,P22);
else
    for ii = 1:Model.NModal
        P11(1,ii) = legendre(X,ii-1);
    end
    P = P11;
end

for ii = 1:Model.NDof 
    Pc{ii,1} = P; 
end

P = blkdiag(Pc{:})+ 1e-26*X;

function y = z1(x)
y = 2*max((sign(-2.0*x+1.0)),0.0)*x;
end

function y = z2(x)
y = max((2.0*(x)-1.0),0.0);
end

function y = w1(x)
if(x>0.5), y=0;
else, y = 1;
end
end

function y = w2(x)
if(x<=0.5), y=0;
else, y = 1;
end
end

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

end
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
