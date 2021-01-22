classdef Model

    properties (Access = public)
        Dof;
        NDof;
        NModal;
        NDisc;
        Table;
        g;
        ge, etae;
        gd;
        tau;
        Texture;
        q;
        dq;
        q0;
        dq0;
        u0;      
        t;
        H;
        Gain;
        Lambda;
        Point;
        Controller;
        Spring;
        Area, Jxx, Jyy, Jzz;
        Xspline;
    end
    
    properties (Access = private)
        Nq;

        E, Nu, Mu;
        xia0;
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
        ActuationSpace;
        AssignProperties;
        
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
    obj.Movie      = false;
    obj.MovieStart = false;
    obj.MovieAxis  = [];
    
    obj.NModal    = 2;
    obj.NDisc     = 2;
    obj.SpaceStep = 11;
    obj.TimeStep  = 1e-1;
    obj.Tdomain   = 15;
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
    obj.xia0    = [0,0,0,1,0,0].';
    obj.Spring  = [1e-9,1];
    
    obj.ActuationSpace = -1;
    obj.AssignProperties = -1;
    
    obj.Point = [1,0,0,0,1,0,0];
    obj.Gain = [1e-2,0.05];
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    cdsoro;
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
    Model = inertia(Model);
    Model = GenerateCosserat(Model);
    Model = GenerateConfigFile(Model);
    
    if isempty(Model.Xspline)
        T = Model.Tdomain;
        X1 = Model.Point(5);
        Model.Xspline = [0, X1; 
                         T, X1 + 1e-6];
    end
end

%--------------------------------------------------------------------- make
function Model = make(Model)

if Model.NDisc > 1
    dir_path = './src/model/tools/solver/build';
else
    dir_path = './src/model/tools/solver/build_continious';
end

tic
CMD = ['cd ',dir_path, ' && make'];
cout('green',['* buidling current dir: ', dir_path, '\n']); pause(1);
system(CMD);
toc    
end

%---------------------------------------------------------------------- set
function Model = csolve(Model)   

if Model.NDisc > 1
    dir_path = './src/model/tools/solver/build_discontinious';
else
    dir_path = './src/model/tools/solver/build_continious';
end

%//////////////////////////////////
writematrix(Model.q0,[dir_path,'/log/state.txt'],'Delimiter','tab');
writematrix(Model.dq0,[dir_path,'/log/momenta.txt'],'Delimiter','tab');
writematrix([zeros(3,1);Model.Gravity(:)],...
    [dir_path,'/log/grav_vector.txt'],'Delimiter','tab');
writematrix(Model.xia0(:),[dir_path,'/log/xi0_vector.txt'],'Delimiter','tab');
writematrix(Model.Xspline,[dir_path,'/log/splineXd.txt'],'Delimiter','tab');

%//////////////////////////////////
out = fullfile(dir_path,'log/tau_vector.log');
fileID = fopen(out,'w');
fprintf(fileID,'%d\n',Model.u0);
fclose(fileID);
%//////////////////////////////////
%out = fullfile(dir_path,'log/splineXd.log');
% fileID = fopen(out,'w');
% fprintf(fileID,'%d\n',Model.Xspline);
% fclose(fileID);
%dlmwrite([dir_path,'/log/splineXd.log'],Model.Xspline);

% solving via c++ executable
tic
CMD = ['cd ',dir_path, ' && ./solver config.m'];
system(CMD);
toc;

% extracting data from log-files
y = load([dir_path ,'/log/state.txt']);
Model.q = y(:,2:end);
Model.t = y(:,1);

y = load([dir_path ,'/log/hamiltonian.txt']);
Model.H = [y(:,2:end), sum(y(:,2:end),2)];

y = load([dir_path ,'/log/endeffector.txt']);
Model.ge = y(:,2:end);

y = load([dir_path ,'/log/endeffector_Vel.txt']);
Model.etae = y(:,2:end);

y = load([dir_path ,'/log/setpoint.txt']);
Model.gd = y(:,2:end);
end

%----------------------------------------------------------------- simulate
function [g, X] = string(Model,k,Quality)
    
    if nargin < 3
        Quality = 100;
    end
    
    N = Model.Nq; 
    Qtmp = Model.q(k,1:N).';
    X = linspace(0,Model.Sdomain,Quality);
    ee = Model.Ba*Model.Phi(Model.Length)*Qtmp;
    
    Model.Length = ee(4);
    
    g0 = [1,0,0,0,0,0,0];
    eta0 = [0,0,0,0,0,0];
    deta0 = [0,0,0,0,0,0];
    
    [~,yf] = ode23(@(t,x) ForwardODE(Model,t,x,Qtmp,...
         Qtmp*0, Qtmp*0),X,[g0 eta0 deta0]);

    g = yf(:,1:7);
end

%----------------------------------------------------------------- simulate
function Model = inertia(Model)
    Model.Area    = pi*Model.Radius^2;
    Model.Jxx     = 0.5*Model.Radius^2;
    Model.Jyy     = 0.25*Model.Radius^2;
    Model.Jzz     = 0.25*Model.Radius^2; 
end

%------------------------------------------------------------ show 3D model
function showModel(Model)
    
figure(101); 
hold all; 
N = Model.Nq;
%Model.Length = Model.Length0;
X = linspace(0,Model.Sdomain,200);

msh = Gmodel('SlenderRod.stl');
%msh = Gmodel('SoftActuatorRedux.stl');
%msh = Gmodel('SoftActuatorPlanarRedux.stl');
%msh = Gmodel('Pneulink.stl'); 
%msh = Gmodel('Pneunet.stl'); 
assignin('base','msh',msh);
mshgr = Gmodel('SoftGripperRedux.stl'); assignin('base','mshgr',mshgr);
%mshgr = Gmodel([]);
%msh = msh.set('Node0',mshgr.Node);
msh = Blender(msh,'Scale',{'x',3});
mshgr = mshgr.set('Node0',mshgr.Node*Model.Sdomain*1e-6);
mshgr = mshgr.set('Node',mshgr.Node*Model.Sdomain*1e-6);
msh = msh.set('Node0',msh.Node*Model.Sdomain);
msh = msh.set('Node',msh.Node*Model.Sdomain);


% set texture
msh.Texture = Model.Texture;
mshgr.Texture = Model.Texture;

msh   = msh.bake();
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
view(-45-90,10);
msh.update();
mshgr.update();
msh.ground(Model.MovieAxis);

if length(Model.t) > 2, FPS = max(round((1/12)/(mean(diff(Model.t)))),1); 
else, FPS = 1; 
end

plotvector([0;0;0],[0.2;0;0],'Color',col(1),'linewidth',2,'MaxHeadSize',0.75);
plotvector([0;0;0],[0;0.2;0],'Color',col(2),'linewidth',2,'MaxHeadSize',0.75);
plotvector([0;0;0],[0;0;0.2],'Color',col(3),'linewidth',2,'MaxHeadSize',0.75);

if ~isempty(Model.xd)
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
        
%         P = Model.Phi;
%         xi_ = zeros(length(X),1);
%         
%         for ii = 1:length(X)
%             xi_(ii) = sum((P(X(ii))*Qtmp).^2); 
%         end
        
        h = plot3(SweepSE3(:,7),SweepSE3(:,6),-SweepSE3(:,5),'linewidth',2,...
            'color',col(2));
%         h = patch(SweepSE3(:,7),SweepSE3(:,6),-SweepSE3(:,5),xi_,'FaceColor',...
%             'none','EdgeColor','interp');

    end
      
    %mshgr.updateNode();
    %msh.updateNode();
    msh.update();
    mshgr.update();
    
    title(['T = ',num2str(Model.t(ii),3)]);
        
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
Model.u0  = zeros(Model.NDof*Model.NModal,1);
end

%---------------------------------------------------------------------- set
function Model = GenerateConfigFile(Model)
   
global FID;

if Model.NDisc > 1
    dir_path = './src/model/tools/solver/build_discontinious';
else
    dir_path = './src/model/tools/solver/build_continious';
end

File = [dir_path,'/config.m'];
if exist(File,'file')
    delete(File);
end

FID = fopen(File,'w');

fprintf(FID,'[options] \n');
if(Model.Controller == 1)
    fprintf(FID,'KINEMATIC_CONTROLLER = 0 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 1 \n');
elseif(Model.Controller == -1)
    fprintf(FID,'KINEMATIC_CONTROLLER = 1 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 0 \n');
else
    fprintf(FID,'KINEMATIC_CONTROLLER = 0 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 0 \n');
end
fprintf(FID,'WRITE_OUTPUT         = 1 \n');
fprintf(FID,['ACTUSPACE =', num2str(Model.ActuationSpace), '\n']);
fprintf(FID,['PROPERTYSET =', num2str(Model.AssignProperties), '\n']);

%fprintf(FID,['GRAV_VECTOR =', num2str(Model.GravityDirection), '\n']);

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
fprintf(FID,['NDOF    = ', num2str(sum(Model.Table)), '\n']);

fprintf(FID,'\n[solver] \n');
fprintf(FID,['TDOMAIN = ', num2str(Model.Tdomain), '\n']);
fprintf(FID,['SPACESTEP = ', num2str(Model.SpaceStep), '\n']);
fprintf(FID,['TIMESTEP  = ', num2str(Model.TimeStep), '\n']);
fprintf(FID,['INTSTEP   = ', '100', '\n']);
fprintf(FID,['ATOL      = ', '1e-2', '\n']);
fprintf(FID,['RTOL      = ', '1e-2', '\n']);
fprintf(FID,['MAX_IMPL  = ', '-1', '\n']);
fprintf(FID,['MAX_ITER  = ', '1', '\n']);
fprintf(FID,['MAX_IK    = ', '1500', '\n']);
fprintf(FID,['SPEEDUP   = ', '80', '\n']);
fprintf(FID,['ADAMPING  = ', '1', '\n']);

fprintf(FID,'\n[physics] \n');
fprintf(FID,['LENGTH   = ', num2str(Model.Sdomain), '\n']);
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
fprintf(FID,['KF1 = ', num2str(Model.Spring(1)), '\n']);
fprintf(FID,['KF2 = ', num2str(Model.Spring(2)), '\n']);
fprintf(FID,['LAMBDA    = ', num2str(Model.Lambda), '\n']);
fprintf(FID,['SPLINEORDER    = ',num2str(1),'\n']);

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
X = single(X)/Model.Sdomain;

if Model.NDisc>1
    P11 = zeros(1,Model.NModal/Model.NDisc);
    P22 = zeros(1,Model.NModal/Model.NDisc);
else
    P11 = zeros(1,Model.NModal);
end


if Model.NDisc>1
    for ii = 1:(Model.NModal/Model.NDisc)
         P11(1,ii) = w1(X)*sign(z1(X)).*legendre(z1(X),ii-1);
         P22(1,ii) = w2(X)*sign(z2(X)).*legendre(z2(X),ii-1);
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
ee = Model.Ba*Model.Phi(t)*q + Model.xia0(:);
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
