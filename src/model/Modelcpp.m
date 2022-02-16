classdef Modelcpp

    properties (Access = public)
        Dof;
        NDof;
        NModal;
        NDisc;
        Table;
        g;
        ge, etae, detae;
        gd;
        tau;
        Texture;
        p,p_;
        q,q_;
        dq;
        q0, q0_;
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
        Kee;
        Xspline;
        Tdomain;
        Sdomain
        Length;
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
%-------------------------------------------------------------- Modelcpp Class
function obj = Modelcpp(Table,varargin) 
    
    obj.Table      = Table;
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
    obj.Gain = [1e-2,0.05,0.02];
    obj.Lambda = [1e-4,1e-4,1e-4];
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    cdsoro;
end
%---------------------------------------------------------------------- get     
function varargout = get(Modelcpp,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Modelcpp.(varargin{ii});
        end
    else
        varargout = Modelcpp.(varargin);
    end
end        
%---------------------------------------------------------------------- set
function Modelcpp = set(Modelcpp,varargin)
    for ii = 1:2:length(varargin)
        Modelcpp.(varargin{ii}) = varargin{ii+1};
    end
end
%---------------------------------------------------------------------- set
function Modelcpp = generate(Modelcpp)
    Modelcpp.Length0 = Modelcpp.Length;
    Modelcpp = inertia(Modelcpp);
    Modelcpp = GenerateCosserat(Modelcpp);
    Modelcpp = GenerateConfigFile(Modelcpp);
    
    Modelcpp = linearStiffness(Modelcpp);
    
    if isempty(Modelcpp.Xspline)
        T = Modelcpp.Tdomain;
        X1 = Modelcpp.Point(5);
        Modelcpp.Xspline = [0, X1; 
                         T, X1 + 1e-6];
    end
end
%--------------------------------------------------------------------- make
function Modelcpp = make(Modelcpp)

%check if cmake present
origLD = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH','');
[status,~] = system("cmake --version");
if  status~= 0
    error('Cmake not installed! Cannot make file')
else
    cout('green','* Cmake compiler founded\n')
end

[status1,~] = system("gcc --version");
[status2,~] = system("g++ --version");
[status3,~] = system("clang --version");

if (status1 ~= 0 && status2 ~= 0 && status3 ~= 0) 
    error('No C++ compiler found! please make sure you have a c-compiler!')
else
    cout('green','* C++ compiler founded\n')
end
    
if Modelcpp.NDisc > 1
    dir_path = './src/Model/tools/solver/build_discontinious';
else
    dir_path = './src/Model/tools/solver/build_continious';
end

tic;
CMD = ['cd ',dir_path, ' && make'];
cout('green',['* building current dir: ', dir_path, '\n']); pause(1);
system(CMD);
toc;

end
%---------------------------------------------------------------------- set
function Modelcpp = csolve(Modelcpp)   

if Modelcpp.NDisc > 1
    dir_path = './src/model/tools/solver/build_discontinious';
else
    dir_path = './src/model/tools/solver/build_continious';
end

%//////////////////////////////////
writeMatrixFile([dir_path,'/log/state.txt'],Modelcpp.q0,'delimiter','\t');
writeMatrixFile([dir_path,'/log/momenta.txt'],Modelcpp.dq0,'delimiter','\t');
writeMatrixFile([dir_path,'/log/estimate.txt'],Modelcpp.q0_,'delimiter','\t');
writeMatrixFile([dir_path,'/log/grav_vector.txt'],...
    [zeros(3,1);Modelcpp.Gravity(:)],'delimiter','\t');
writeMatrixFile([dir_path,'/log/xi0_vector.txt'],Modelcpp.xia0(:),'delimiter','\t');
writeMatrixFile([dir_path,'/log/splineXd.txt'],Modelcpp.Xspline,'delimiter','\t');

%//////////////////////////////////
out = fullfile(dir_path,'log/tau_vector.log');
fileID = fopen(out,'w');
fprintf(fileID,'%d\n',Modelcpp.u0);
fclose(fileID);
%//////////////////////////////////

% solving via c++ executable
tic
CMD = ['cd ',dir_path, ' && ./solver config.m'];
system(CMD);
toc;

% extracting data from log-files
y = load([dir_path ,'/log/state.txt']);
Modelcpp.q = y(:,2:end);
Modelcpp.t = y(:,1);

% extracting data from log-files
y = load([dir_path ,'/log/momenta.txt']);
Modelcpp.p = y(:,2:end);

y = load([dir_path ,'/log/estimate.txt']);
Modelcpp.q_ = y(:,2:end);

y = load([dir_path ,'/log/hamiltonian.txt']);
Modelcpp.H = [y(:,2:end), sum(y(:,2:end),2)];

y = load([dir_path ,'/log/endeffector.txt']);
Modelcpp.ge = y(:,2:end);

y = load([dir_path ,'/log/endeffector_Vel.txt']);
Modelcpp.etae = y(:,2:end);

y = load([dir_path ,'/log/setpoint.txt']);
Modelcpp.gd = y(:,2:end);

y = load([dir_path ,'/log/endeffector_Acc.txt']);
Modelcpp.detae = y(:,2:end);
end
%-------------------------------------------------- compute Cosserat string
function [g, X] = string(Modelcpp,z,varargin)
    
%     if nargin < 3
      Quality = 200;
%     end
%     
%     if nargin < 2
%         Qtmp = Modelcpp.q0;
%     else
%         if k <= 0
%             Qtmp = Modelcpp.q0;
%         else
%             Qtmp = Modelcpp.q(k,1:Modelcpp.Nq).';
%         end
%     end
%     
%     if isempty(varargin)
%         %s = linspace(0,Modelcpp.Sdomain,Quality);
%     else
%         Qtmp = varargin{1};
%     end
%     
    Qtmp = z(:);
    X = linspace(0,Modelcpp.Sdomain,Quality);
    ee = Modelcpp.Ba*Modelcpp.Phi(Modelcpp.Length)*Qtmp;
    
    Modelcpp.Length = ee(4);
    
    g0 = [1,0,0,0,0,0,0];
    eta0 = [0,0,0,0,0,0];
    deta0 = [0,0,0,0,0,0];
    
    [~,yf] = ode23(@(t,x) ForwardODE(Modelcpp,t,x,Qtmp,...
         Qtmp*0, Qtmp*0),X,[g0 eta0 deta0]);

    g = yf(:,1:7);

end
%------------------------------------------- compute Cosserat strain fields
function [K,s] = strainField(Modelcpp,k,Quality,varargin)
    
if nargin < 3
    Quality = 100;
end

if nargin < 2
    Qtmp = Modelcpp.q0;
else
    if k <= 0
        Qtmp = Modelcpp.q0;
    else
        Qtmp = Modelcpp.q(k,1:Modelcpp.Nq).';
    end
end
    
if isempty(varargin)
    s = linspace(0,Modelcpp.Sdomain,Quality);
else
    s = varargin{1};
end

K = zeros(length(s),Modelcpp.NDof);

for jj = 1:length(s)
    K(jj,:) = transpose((Modelcpp.Phi(s(jj))*Qtmp));
end
    
end

%--------------------------------------------------- update animation frame
function Modelcpp = updateFrame(Modelcpp)
   
    if Modelcpp.Movie
        background(gitpage);
        if Modelcpp.MovieStart == false
            Modelcpp.MovieStart = true;
            MovieMaker(Modelcpp,'Modelcpp','Start');
        else
            MovieMaker(Modelcpp,'Modelcpp','');
        end
    end
    
end
%-------------------------------------------------- compute inertia of disk 
function Modelcpp = inertia(Modelcpp)
    Modelcpp.Area    = pi*Modelcpp.Radius^2;
    Modelcpp.Jxx     = 0.5*Modelcpp.Radius^2;
    Modelcpp.Jyy     = 0.25*Modelcpp.Radius^2;
    Modelcpp.Jzz     = 0.25*Modelcpp.Radius^2; 
end
%------------------------------------------------- initial linear stiffness 
function Modelcpp = linearStiffness(Modelcpp)
    E0 = Modelcpp.E;
    G0 = E0/(2*(1+Modelcpp.Nu));
    A = Modelcpp.Area;
    J11 = Modelcpp.Jxx;
    J22 = Modelcpp.Jyy;
    J33 = Modelcpp.Jzz;
    
    S = diag([G0*J11,E0*J22,E0*J33,E0*A,G0*A,G0*A]);
    s = linspace(0,Modelcpp.Sdomain,150);
    ds = mean(diff(s));
    
    Modelcpp.Kee = zeros(Modelcpp.NDof*Modelcpp.NModal,Modelcpp.NDof*Modelcpp.NModal);
    for ii = 1:Modelcpp.SpaceStep
        Modelcpp.Kee = Modelcpp.Kee + ds*(Modelcpp.Ba*Modelcpp.Phi(s(ii))).'*S*Modelcpp.Ba*Modelcpp.Phi(s(ii));
    end
end
%------------------------------------------------------------ show 3D Modelcpp
function showModelcpp(Modelcpp)
    
figure(101); 
hold all; 
N = Modelcpp.Nq;
%Modelcpp.Length = Modelcpp.Length0;
X = linspace(0,Modelcpp.Sdomain,200);

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
mshgr = mshgr.set('Node0',mshgr.Node*Modelcpp.Sdomain*1e-6);
mshgr = mshgr.set('Node',mshgr.Node*Modelcpp.Sdomain*1e-6);
msh = msh.set('Node0',msh.Node*Modelcpp.Sdomain);
msh = msh.set('Node',msh.Node*Modelcpp.Sdomain);


% set texture
msh.Texture = Modelcpp.Texture;
mshgr.Texture = Modelcpp.Texture;

msh   = msh.bake();
mshgr = mshgr.bake();
msh   = msh.render();
mshgr = mshgr.render();

axis equal; 

LinkID = knnsearch(X(:),msh.Node(:,3));
if isempty(Modelcpp.MovieAxis)
    Modelcpp.MovieAxis = boxhull(msh.Node); 
end
axis(Modelcpp.MovieAxis);
drawnow;
view(-45-90,10);
msh.update();
mshgr.update();
msh.ground(Modelcpp.MovieAxis);

if length(Modelcpp.t) > 2, FPS = max(round((1/12)/(mean(diff(Modelcpp.t)))),1); 
else, FPS = 1; 
end

plotvector([0;0;0],[0.2;0;0],'Color',col(1),'linewidth',2,'MaxHeadSize',0.75);
plotvector([0;0;0],[0;0.2;0],'Color',col(2),'linewidth',2,'MaxHeadSize',0.75);
plotvector([0;0;0],[0;0;0.2],'Color',col(3),'linewidth',2,'MaxHeadSize',0.75);

if ~isempty(Modelcpp.xd)
p = Modelcpp.xd(end,5:7);
R = quat2rot(Modelcpp.xd(end,1:4));
fr = [p(1,3);p(1,2);-p(1,1)];
% % 
f0 = plotvector(fr,R*[0.2;0;0],'Color',col(1),'linewidth',2,'MaxHeadSize',0.75);
f1 = plotvector(fr,R*[0;0.2;0],'Color',col(2),'linewidth',2,'MaxHeadSize',0.75);
f2 = plotvector(fr,R*[0;0;0.2],'Color',col(3),'linewidth',2,'MaxHeadSize',0.75);
end

vline = [];

for ii = [1:FPS:length(Modelcpp.t),length(Modelcpp.t)]

    msh.reset();
    mshgr.reset();

    Qtmp = Modelcpp.q(ii,1:N).';
    ee = Modelcpp.Ba*Modelcpp.Phi(Modelcpp.Length)*Qtmp;
    
    Modelcpp.Length = ee(4);
    
    g0 = [1,0,0,0,0,0,0];
    eta0 = [0,0,0,0,0,0];
    deta0 = [0,0,0,0,0,0];
    
    [~,yf] = ode23(@(t,x) ForwardODE(Modelcpp,t,x,Qtmp,...
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
       
    if ~isempty(Modelcpp.xd)
       p = Modelcpp.xd(:,4:6);
       h0 = plot3(p(:,3),p(:,2),-p(:,1),'k--');
    end
        
    else
        
        delete(h);
         if ~isempty(Modelcpp.xd)
        delete(hm);
        delete(h0);
        end
       if size(vline,1) > 4, delete(hvm); end
       
        
        if ~isempty(Modelcpp.xd)
        delete(h0);
        delete(f0);
        delete(f1);
        delete(f2);
        %delete(h1);
        %delete(h2);
        p = Modelcpp.xd(ii,5:7);
        R = quat2rot(Modelcpp.xd(ii,1:4));
        fr = [p(3);p(2);-p(1)];

         f0 = plotvector(fr,R*[0.2;0;0],'Color',col(1),'linewidth',2,'MaxHeadSize',0.75);
         f1 = plotvector(fr,R*[0;0.2;0],'Color',col(2),'linewidth',2,'MaxHeadSize',0.75);
         f2 = plotvector(fr,R*[0;0;0.2],'Color',col(3),'linewidth',2,'MaxHeadSize',0.75);
        end
        
%         P = Modelcpp.Phi;
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
    
    title(['T = ',num2str(Modelcpp.t(ii),3)]);
        
    if ~isempty(Modelcpp.MovieAxis), axis(Modelcpp.MovieAxis); end
    
    if Modelcpp.Movie
        if Modelcpp.MovieStart == false
            Modelcpp.MovieStart = true;
            MovieMaker(Modelcpp,'mdl','Start');
        else
            MovieMaker(Modelcpp);
        end
    end
end
    
end

end
%--------------------------------------------------------------------------
methods (Access = private)

%---------------------------------------------------------------------- set
function Modelcpp = GenerateCosserat(Modelcpp)
I6 = eye(6);
set = 1:6;
xa = set(logical(Modelcpp.Table));
xc = setdiff(set,xa);

Modelcpp.Ba   = I6(:,xa);
Modelcpp.Bc   = I6(:,xc);
Modelcpp.NDof = size(Modelcpp.Ba,2);
Modelcpp.Phi  = @(x) ShapeFunction(Modelcpp,x);
Modelcpp.Nq   = Modelcpp.NDof*Modelcpp.NModal;
Modelcpp.q0   = zeros(Modelcpp.NDof*Modelcpp.NModal,1);
Modelcpp.q0_  = zeros(Modelcpp.NDof*Modelcpp.NModal,1);
Modelcpp.dq0  = zeros(Modelcpp.NDof*Modelcpp.NModal,1);
Modelcpp.u0   = zeros(Modelcpp.NDof*Modelcpp.NModal,1);
end

%---------------------------------------------------------------------- set
function Modelcpp = GenerateConfigFile(Modelcpp)
   
global FID;

if Modelcpp.NDisc > 1
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
if(Modelcpp.Controller == 1)
    fprintf(FID,'KINEMATIC_CONTROLLER = 0 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 1 \n');
elseif(Modelcpp.Controller == -1)
    fprintf(FID,'KINEMATIC_CONTROLLER = 1 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 0 \n');
else
    fprintf(FID,'KINEMATIC_CONTROLLER = 0 \n');
    fprintf(FID,'ENERGY_CONTROLLER    = 0 \n');
end
fprintf(FID,'WRITE_OUTPUT         = 1 \n');
fprintf(FID,['ACTUSPACE =', num2str(Modelcpp.ActuationSpace), '\n']);
fprintf(FID,['PROPERTYSET =', num2str(Modelcpp.AssignProperties), '\n']);

%fprintf(FID,['GRAV_VECTOR =', num2str(Modelcpp.GravityDirection), '\n']);

if(Modelcpp.NDisc > 1)
    fprintf(FID,'DISCONTINIOUS        = 1 \n');
else
    fprintf(FID,'DISCONTINIOUS        = 0 \n');
end

fprintf(FID,'\n[cosserat] \n');
fprintf(FID,['K1 = ', num2str(Modelcpp.Table(1)), '\n']);
fprintf(FID,['K2 = ', num2str(Modelcpp.Table(2)), '\n']);
fprintf(FID,['K3 = ', num2str(Modelcpp.Table(3)), '\n']);
fprintf(FID,['E1 = ', num2str(Modelcpp.Table(4)), '\n']);
fprintf(FID,['E2 = ', num2str(Modelcpp.Table(5)), '\n']);
fprintf(FID,['E3 = ', num2str(Modelcpp.Table(6)), '\n']);

fprintf(FID,'\n[model] \n');
fprintf(FID,['NMODE   = ', num2str(Modelcpp.NModal), '\n']);
fprintf(FID,['NDISC   = ', num2str(Modelcpp.NDisc), '\n']);
fprintf(FID,['NDOF    = ', num2str(sum(Modelcpp.Table)), '\n']);

fprintf(FID,'\n[solver] \n');
fprintf(FID,['TDOMAIN = ', num2str(Modelcpp.Tdomain), '\n']);
fprintf(FID,['SPACESTEP = ', num2str(Modelcpp.SpaceStep), '\n']);
fprintf(FID,['TIMESTEP  = ', num2str(Modelcpp.TimeStep), '\n']);
fprintf(FID,['INTSTEP   = ', '100', '\n']);
fprintf(FID,['ATOL      = ', '1e-2', '\n']);
fprintf(FID,['RTOL      = ', '1e-2', '\n']);
fprintf(FID,['MAX_IMPL  = ', '-1', '\n']);
fprintf(FID,['MAX_ITER  = ', '1', '\n']);
fprintf(FID,['MAX_IK    = ', '1500', '\n']);
fprintf(FID,['SPEEDUP   = ', '80', '\n']);
fprintf(FID,['ADAMPING  = ', '1', '\n']);

fprintf(FID,'\n[physics] \n');
fprintf(FID,['LENGTH   = ', num2str(Modelcpp.Sdomain), '\n']);
fprintf(FID,['RHO      = ', num2str(Modelcpp.Density), '\n']);
fprintf(FID,['EMOD     = ', num2str(Modelcpp.E), '\n']);
fprintf(FID,['NU       = ', num2str(Modelcpp.Nu), '\n']);
fprintf(FID,['MU       = ', num2str(Modelcpp.Mu), '\n']);
fprintf(FID,['PRS_AREA = ', '1e-5', '\n']);
fprintf(FID,['GRAVITY  = ', num2str(Modelcpp.Gravity), '\n']);
fprintf(FID,['RADIUS   = ', '0.01', '\n']);
fprintf(FID,['AREA     = ', num2str(Modelcpp.Area), '\n']); %'3.0000e-04'
fprintf(FID,['J_XX     = ', num2str(Modelcpp.Jxx), '\n']); % '2.7250e-10'
fprintf(FID,['J_YY     = ', num2str(Modelcpp.Jyy), '\n']); % '2.2500e-11'
fprintf(FID,['J_ZZ     = ', num2str(Modelcpp.Jzz), '\n']); %'2.5000e-10'

fprintf(FID,'\n[control] \n');
fprintf(FID,['KP = ', num2str(Modelcpp.Gain(1)), '\n']);
fprintf(FID,['KD = ', num2str(Modelcpp.Gain(2)), '\n']);
fprintf(FID,['LK = ', num2str(Modelcpp.Gain(3)), '\n']);
fprintf(FID,['KF1 = ', num2str(Modelcpp.Spring(1)), '\n']);
fprintf(FID,['KF2 = ', num2str(Modelcpp.Spring(2)), '\n']);
fprintf(FID,['LAMBDA    = ', num2str(Modelcpp.Lambda(1)), '\n']);
fprintf(FID,['LAMBDAK   = ', num2str(Modelcpp.Lambda(2)), '\n']);
fprintf(FID,['SPLINEORDER    = ',num2str(1),'\n']);

fprintf(FID,'\n[setpoint] \n');
fprintf(FID,['Q1d = ', num2str(Modelcpp.Point(1)), '\n']);
fprintf(FID,['Q2d = ', num2str(Modelcpp.Point(2)), '\n']);
fprintf(FID,['Q3d = ', num2str(Modelcpp.Point(3)), '\n']);
fprintf(FID,['Q4d = ', num2str(Modelcpp.Point(4)), '\n']);
fprintf(FID,['Xd =  ', num2str(Modelcpp.Point(5)), '\n']);
fprintf(FID,['Yd =  ', num2str(Modelcpp.Point(6)), '\n']);
fprintf(FID,['Zd =  ', num2str(Modelcpp.Point(7)), '\n']);

fclose('all');
    
end

%---------------------------------------------------------------------- set
function P = ShapeFunction(Modelcpp,X)

Pc = cell(Modelcpp.NDof,1);
X = single(X)/Modelcpp.Sdomain;

if Modelcpp.NDisc>1
    P11 = zeros(1,Modelcpp.NModal/Modelcpp.NDisc);
    P22 = zeros(1,Modelcpp.NModal/Modelcpp.NDisc);
else
    P11 = zeros(1,Modelcpp.NModal);
end


if Modelcpp.NDisc>1
    for ii = 1:(Modelcpp.NModal/Modelcpp.NDisc)
         P11(1,ii) = w1(X)*sign(z1(X)).*legendre(z1(X),ii-1);
         P22(1,ii) = w2(X)*sign(z2(X)).*legendre(z2(X),ii-1);
         %P11(1,ii) = w1(X)*chebyshev(z1(X),ii-1);
         %P22(1,ii) = w2(X)*chebyshev(z2(X),ii-1);

    end
    P = horzcat(P11,P22);
else
    for ii = 1:Modelcpp.NModal
        P11(1,ii) = legendre(X,ii-1);
    end
    P = P11;
end

for ii = 1:Modelcpp.NDof 
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
function dg = ForwardODE(Modelcpp,t,g,q,dq,ddq)
ee = Modelcpp.Ba*Modelcpp.Phi(t)*q + Modelcpp.xia0(:);
dee = Modelcpp.Ba*Modelcpp.Phi(t)*dq;
ddee = Modelcpp.Ba*Modelcpp.Phi(t)*ddq;

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
function MovieMaker(Modelcpp,Name,Request)
if nargin < 2, Request = ''; end

if Modelcpp.Movie
    switch(Request)
        case('Start')
            filename = string([Name,'_', char(datetime(now,...
                'ConvertFrom','datenum')),'.gif']);
            
            filename = erase(filename,[":"," "]);
            background(metropolis);
            if ~isempty(Modelcpp.MovieAxis), axis(Modelcpp.MovieAxis); end
            drawnow;
            gif(char(filename),'frame',gcf,'nodither');
        otherwise
            background(metropolis);
            if ~isempty(Modelcpp.MovieAxis), axis(Modelcpp.MovieAxis); end
            drawnow;
            gif;
    end
end
end
