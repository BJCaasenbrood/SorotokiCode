classdef Shapes

    properties (Access = public)
        NDim;       % State dimensions
        NJoint;     % Joint dimension (i.e, NDim/2)
        NDof;       % Shape joint DoF (max 6 -- geometric freedom)
        NInput;     % Number of inputs (=dim(q) by default)
        NModal;     % Number of modes
        NNode;      % Number of nodes (spatial discretization)
        
        Length;     % Intrinic length
        Theta;      % Shape Basis Matrix functional
        Xi0;        % Intrinsic strain fucntional
        g0;         % Initial orientation frame
        gL;         % End-efector frame
        X0;         % Initial conditions x0:=(q0,dq0)
              
        Texture;    % Texture for render
        Material;   % Material model    
        Gravity;    % Gravity vector
        Sdf;        % Cross-sectional SDF
        Contact;    % Contact model
        Muscle;     % Muscle function
        
        Fem;        % a-priori fem solution
        Sigma;      % spatial discretization

        InputMap;   % (Nonlinear) input mapping G(q)
        Log;        % Log file

        Node;       % 
        Node0;     
        
        TubeRadiusA;
        TubeRadiusB;
        TubeRadiusAlpha;
        TubeRamp;
        
        TubeMatCap;
        ContactDistance;
        
        Kp;     % IK proportional gain
        Kd;     % IK differential gain
    end
    
    properties (Access = private)
        Table;
        Gmodel;
        xia0;
        ds;
        Ba;
        
        Att;
        Jtt;
        Mtt;
        Ktt;
        Dtt;
        
        Center;
        Rotation;
        Gamma;
        Kappa;
                
        POD;
        PODR;
        PODQ;
        
        PODEnergy;
        
        Filter;
        FilterRadius;
        Quality;
        
        ThetaEval;
        Xi0Eval;
        MuscleEval;
        VolumetricContact;
        MuscleLines;
        RenderMuscles = true;
        
        Umin = 0;
        Umax = 1;
    end
   
%--------------------------------------------------------------------------
methods  
%----------------------------------------------- MODAL SHAPE RECONSTRUCTION
function obj = Shapes(Input,NModal,varargin) 
    
    obj.NModal = NModal;
    obj.Table  = double(NModal > 0); 
    obj.NDof   = sum(obj.Table);
    obj.xia0   = [0,0,0,0,0,1].';
    
    if isa(Input,'Fem')
        obj.Fem    = Input;
        obj.NNode  = 30;        
        gvec = Input.get('Gravity');
        if ~isempty(gvec)
            if numel(gvec) == 2
                obj.Gravity = [0; gvec(:)];
            else
                obj.Gravity = gvec(:);
            end
        else
            obj.Gravity = zeros(3,1);
        end
        
    elseif isa(Input,'double')
        obj.NNode = size(Input,1);
        obj.PODQ = Input;
        obj.PODR = Input;
        
        obj.Length    = 100;
        obj.Sigma = linspace(0,1,obj.NNode); 
        obj.ds    = obj.Length/(obj.NNode);

        obj.PODEnergy{2} = ones(size(Input,2),1);
        obj.PODEnergy{1} = ones(size(Input,2),1);
       
        obj.Gravity = zeros(3,1);
    end

    % cross-section SDF
    obj.Sdf          = sCircle(5);
    obj.Center       = zeros(3,1);
    obj.Material     = NeoHookeanMaterial(5,0.33);
    obj.Texture      = bluebase;
    
    obj.FilterRadius      = 10;
    obj.VolumetricContact = true;
    obj.ContactDistance   = 1e-3;

    obj.g0 = SE3(eye(3),zeros(3,1));
    obj.gL = [];
    obj.Log = struct;
    
    obj.Log.FK  = [];
    obj.Log.EL  = [];
    obj.Log.Shp = [];

    obj.Kp = 1e-4;
    obj.Kd = 1e3;
    obj.TubeRadiusA     = 5;
    obj.TubeRadiusB     = 5;
    obj.TubeRadiusAlpha = 0;  
    obj.TubeRamp        = 0;  
    obj.Umin            = -1;
    obj.Umax            = 1;
           
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
   
    if ~isempty(obj.Fem)
    if ~isempty(obj.Fem.get('Output'))
        out = obj.Fem.get('Output');
        Nd0 = obj.Fem.get('Node0');
        
        BdBox = boxhull(Nd0(out(:,1),:));
        obj   = reference(obj,[BdBox(1),BdBox(3)],...
            [BdBox(2),BdBox(2)]);
    end
    end

    obj.NJoint = sum(obj.NModal);
    obj.NDim  = 2*obj.NJoint;

    if isempty(obj.X0)
        obj.X0 = zeros(obj.NDim,1);
    end
    
    if isempty(obj.InputMap)
        obj.InputMap = @(x) eye(obj.NJoint);
        obj.NInput = obj.NJoint;
    end
    
    obj = rebuild(obj);
    
end
%---------------------------------------------------------------------- get     
function varargout = get(Shapes,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Shapes.(varargin{ii});
        end
    else
        varargout = Shapes.(varargin);
    end
end       
%---------------------------------------------------------------------- set
function Shapes = set(Shapes,varargin)
    for ii = 1:2:length(varargin)
        Shapes.(varargin{ii}) = varargin{ii+1};
    end
end
%----------------------------------------------------------- set base SE(3)
function Shapes = setBase(Shapes,varargin)
if numel(varargin{1}) == 9 || numel(varargin{1}) == 3
    Shapes.g0 = SE3(varargin{1});
else
    Shapes.g0 = varargin{1};
end     
end
%--------------------------------------------------------------- set radius
function Shapes = setRadius(Shapes,varargin)
   if numel(varargin) == 1 && numel(varargin{1}) == 1
       Shapes.TubeRadiusA = varargin{1};
       Shapes.TubeRadiusB = varargin{1};
   elseif numel(varargin) == 1 && numel(varargin{1}) == 2
      R = varargin{1};
      Shapes.TubeRadiusA = R(1);
      Shapes.TubeRadiusB = R(2); 
   elseif numel(varargin) == 1 && numel(varargin{1}) == 3
      R = varargin{1};
      Shapes.TubeRadiusA = R(1);
      Shapes.TubeRadiusB = R(2);  
      Shapes.TubeRadiusAlpha = R(3);
   end
   
   Shapes.Sdf = sCircle(Shapes.TubeRadiusA);
   Shapes = rebuild(Shapes);
end
%----------------------------------------------------------------- set ramp
function Shapes = setRamp(Shapes,x)
    y = clamp(x,0,1);  
    Shapes.TubeRamp = y;
end
%-------------------------------------------------------------- set gravity 
function Shapes = setInputMap(Shapes,varargin)
Shapes.InputMap = varargin{1};
G0 = Shapes.InputMap(Shapes);
Shapes.NInput = size(G0,2);
end
%-------------------------------------------------------------- set gravity 
function Shapes = setMaterial(Shapes,varargin)
   Shapes.Material = varargin{1};
end
%-------------------------------------------------------------- set gravity 
function Shapes = setTexture(Shapes,varargin)
   Shapes.Texture = varargin{1};
end
%-------------------------------------------------------------- set gravity 
function Shapes = addGravity(Shapes,varargin)
if isempty(varargin)
     vec = -9810*Shapes.xia0(4:6);
     varargin{1} = vec(:);
end
Shapes.Gravity = varargin{1};    
end
%-------------------------------------------------------------- set gravity 
function Shapes = addContact(Shapes,varargin)
Shapes.Contact = varargin{1};    
end
%--------------------------------------------------------------- add muscle
function Shapes = addMuscle(Shapes,varargin)
    Shapes.Muscle{end+1} = varargin{1};
end
%---------------------------------------------------------------- show mesh
function Shapes = show(Shapes,varargin)
    
if nargin<2, Request = -1; 
else, Request = varargin{1}; 
end

figure(101);

switch(Request)
    case('Base'),  Request = 'POD';
    case('POD'),   Request = 'POD';
    otherwise,     Request = 'POD';
end

if strcmp(Request,'POD') 
    
    X = Shapes.Sigma*Shapes.Length;
    Y = [];
    
    for ii = 1:numel(X)
       S = diag(Shapes.Theta(X(ii)));
       Y = [Y; S(:).'];
    end
    
    figure(101); clf;
    
    M = Shapes.NModal; 
    Modes = M(M~=0);
    List  = 1:6;
    List  = List(M~=0);
    
    kk = 1;
    for ii = 1:numel(Modes)
       
        subplot(1,Shapes.NDof,ii);

        for jj = 1:Modes(ii)
            shade(X,Y(:,kk),'LineW',2,...
                'Color',col(jj)); hold on;
            kk = kk + 1;
        end
        
        xlim([min(X) max(X)]);
        ax = gca;
        ax.LineWidth  = 1.5;
        ax.XTick      = [min(X) max(X)];
        ax.XTickLabel = {'0','L'};
        axis square; grid on; box on;
        
        switch(List(ii))
            case 1, title(['$\kappa_{xx}$'],'interp','latex','Fontsize',24);
            case 2, title(['$\kappa_{xz}$'],'interp','latex','Fontsize',24);
            case 3, title(['$\kappa_{yz}$'],'interp','latex','Fontsize',24);
            case 4, title(['$\epsilon_{xx}$'],'interp','latex','Fontsize',24);
            case 5, title(['$\epsilon_{xz}$'],'interp','latex','Fontsize',24);
            case 6, title(['$\epsilon_{yz}$'],'interp','latex','Fontsize',24);
        end

    end
    
end

end
%---------------------------------------------------- render muscle lines
function h = showMuscle(Shapes,varargin)
    
% if nargin < 2
%     q = ones(Shapes.NJoint,1)*1e-12;
% else
%     q = varargin{1};
% end

% [G] = string(Shapes,q);
G = Shapes.Log.FK.g;
M = Shapes.MuscleEval;

hold on;
for ii = 1:size(M,3)
    P = M(1:2:end,1:3,ii);
    for jj = 1:Shapes.NNode
        P_ = G(1:3,1:3,jj)*(P(jj,:).') + G(1:3,4,jj);
        PP(jj,:) = P_.';
    end
    h{ii} = fplot(PP);
end

end
%------------------------------------------------------------ render string
function Shapes = render(Shapes,varargin)

    if nargin < 2
       q = ones(Shapes.NJoint,1)*1e-12; 
    else
       q = varargin{1};
    end
    
    p0 = Shapes.g0(1:3,4);
    [g, J] = string(Shapes,q);
    
    Shapes.Node = [p0.'; reshape(g(1:3,4,1:end),3,[]).'];
    Shapes.Log.FK.g = g;
    Shapes.Log.FK.J = J;
    Shapes.Log.FK.Node = Shapes.Node;
    Shapes.Log.FK.L    = sum(sqrt(sum(diff(Shapes.Log.FK.Node).^2,2))) ...
                         + Shapes.Length/Shapes.NNode;

    Stretch = Shapes.Log.FK.L/Shapes.Length;

    if isempty(Shapes.Gmodel) 
        obj = Gmodel(Shapes,varargin{2:end});
        obj = obj.set('Texture',Shapes.Texture);
        obj = obj.bake;
        obj = obj.render();
        Shapes.Gmodel = obj;
    else
        [x,y,z] = rtubeplot(Shapes.Node.',...
            Shapes.TubeRadiusA,...
            Shapes.TubeRadiusB,...
            Shapes.TubeRadiusAlpha,...
            16,1e-6,...
            Shapes.TubeRamp*Stretch^3);

        fv = surf2patch(x,y,z,'triangles');
        v = fv.vertices;
        Shapes.Gmodel.Node = v;
        Shapes.Gmodel.update();
    end
    
    hold on;
    if Shapes.RenderMuscles && ~isempty(Shapes.Muscle) 
        M = Shapes.MuscleEval;
        P = zeros(Shapes.NNode,3,numel(Shapes.Muscle));
        
        for ii = 1:size(M,3)
                        
            Nds = M(1:2:end,1:3,ii);
            for jj = 1:Shapes.NNode
                P(jj,:,ii) = (g(1:3,1:3,jj)*(Nds(jj,:).') ...
                    + g(1:3,4,jj)).';
            end
        end
        
        Ns   = Shapes.NNode;
        Cmap = polarmap(Ns);
        Xmap = linspace(0,1,Ns);
 
        if isempty(Shapes.MuscleLines)
            Shapes.MuscleLines = fplot(P);
        else
            h = Shapes.MuscleLines;
            for ii = 1:numel(h)
               h(ii).XData = P(:,1,ii);
               h(ii).YData = P(:,2,ii);
               h(ii).ZData = P(:,3,ii);
               
               %v = invlerp(U(ii),Shapes.Umin,Shapes.Umax);
               %h(ii).Color = interp1(Xmap.',Cmap,clamp(v,0,1));
            end
        end
        
    end
    
end
%------------------------------------------------------------ set reference
function Shapes = reference(Shapes,varargin)

if numel(varargin) == 2
    % two points curve
    X0 = varargin{1}; XL = varargin{2};
    
    % build Cosserat curve in reference
    Shapes.Node0 = [linspace(X0(1),XL(1),Shapes.NNode).', ...
        linspace(X0(2),XL(2),Shapes.NNode).'];
    
    %B = boxhull(Shapes.Node0);
    Shapes.Length = sqrt((XL(1)-X0(1))^2 + (XL(2)-X0(2))^2);
elseif numel(varargin) == 1 && size(varargin{1},2) == 2
    Shapes.Node0 = varargin{1};
    Shapes.NNode = length(Shapes.Node0);
    Shapes.Length = sum(sqrt(sum(diff(Shapes.Node0).^2,2)),1);
end

% build discretization of curve domain
Shapes.Sigma = linspace(0,Shapes.Length,Shapes.NNode);    
Shapes.ds    = Shapes.Length/(Shapes.NNode);

% finds associated nodes from Fem mesh.
if isempty(Shapes.Fem)
    Shapes = GenerateRadialFilter(Shapes);
end

Shapes = rebuild(Shapes);

end   
%------------------------------------------------------------ set reference
function Shapes = rebuild(Shapes,varargin)
    
for ii = 1:2:length(varargin)
    Shapes.(varargin{ii}) = varargin{ii+1};
end

set = 1:6;
I6  = eye(6);
Xa  = [];

for ii = 1:6
    for jj = 1:Shapes.NModal(ii)
        Xa = [Xa,set(ii)];
    end
end

Shapes.NJoint  = sum(Shapes.NModal);
Shapes.NDim    = 2*Shapes.NJoint;
Shapes.Ba      = I6(:,Xa);
Shapes.Sigma   = linspace(0,1,Shapes.NNode);
Shapes.ds      = Shapes.Length/(Shapes.NNode);

Shapes = BuildInertia(Shapes);
JJ     = Shapes.Mtt/Shapes.Material.Rho;

% linear approximation of the stiffness
[~,KK] = Shapes.Material.PiollaStress(eye(3));
KK     = 4.15*diag(voightextraction(KK));

% E  = Shapes.Material.getModulus();
% Nu = Shapes.Material.Nu;
% G  = E/(2*(Nu+1));
% KK = pi*diag([G,E,E,E,G,G]);

Shapes.Ktt  = diag(diag(JJ))*KK;
%Shapes.Ktt(4:6,4:6) = 0.001*Shapes.Ktt(4:6,4:6);
Shapes.Dtt  = Shapes.Material.Zeta*Shapes.Ktt;
  
if ~isempty(Shapes.Node0)
    Shapes = GenerateRadialFilter(Shapes);
end

if ~isempty(Shapes.PODR) || ~isempty(Shapes.PODQ)
    
    if length(Shapes.PODQ) ~= Shapes.NNode
        X = linspace(0,Shapes.Length,length(Shapes.PODQ));
        Shapes.PODR = interp1(X,Shapes.PODR,Shapes.Sigma*Shapes.Length);
        Shapes.PODQ = interp1(X,Shapes.PODQ,Shapes.Sigma*Shapes.Length);
    end
    
    % ensure orthonormality
    Shapes.PODR = gsogpoly(Shapes.PODR,Shapes.Sigma);
    Shapes.PODQ = gsogpoly(Shapes.PODQ,Shapes.Sigma);
    
    k = 1;
    Shapes.POD = [];
    for ii = 1:numel(Shapes.NModal)
        for jj = 1:Shapes.NModal(ii)
            if ii == 1
                Shapes.POD(:,k) = Shapes.PODR(:,jj);
            else
                Shapes.POD(:,k) = Shapes.PODQ(:,jj);
            end
            k = k+1;
        end
    end
    
    % rebuild shape-function matrix
    Shapes.Theta = @(x) ShapeFunction(Shapes,x);
end

if ~isa(Shapes.Xi0,'function_handle')
    Shapes.Xi0 = @(x) IntrinsicFunction(Shapes,x);   
end

if ~isempty(Shapes.Theta) 
    
% precompute Theta matrix
FncT = @(x) Shapes.Theta(x);
FncX = @(x) Shapes.Xi0(x);
s   = sort([Shapes.Sigma*Shapes.Length,...
            Shapes.Sigma*Shapes.Length+(2/3)*Shapes.ds]);

[nx,ny]          = size(FncT(0));
Shapes.ThetaEval = zeros(nx,ny,numel(s));
Shapes.Xi0Eval   = zeros(6,1,numel(s));

for ii = 1:numel(s)
    Shapes.ThetaEval(:,:,ii) = FncT(s(ii));
    Shapes.Xi0Eval(:,1,ii)   = FncX(s(ii));
end
end

if ~isempty(Shapes.Muscle)
    M = numel(Shapes.Muscle);
    N = numel(s);
    h = mean(diff(s/Shapes.Length));
    
    Shapes.MuscleEval = zeros(N,6,M);
    [P, P0] = computeMuscleGroups(Shapes,s/Shapes.Length);
    
    for ii = 1:M
        
       % compute derivative along the tendon 
       [~,dP0] = gradient(P0(:,:,ii),h); 

       Shapes.MuscleEval(:,1:3,ii) = P0(:,:,ii);
       Shapes.MuscleEval(:,4:6,ii) = dP0;
    end
    
    Shapes.NInput = M;
else
   N = numel(s);
   Shapes.MuscleEval = zeros(N,6,1);
end

end 
%--------------------------------------------------------------------------
function Shapes = reconstruct(Shapes,varargin)
    
fem = Shapes.Fem;
t   = fem.Log.t;

set = 1:6;
I6  = eye(6);
xa  = set(logical(Shapes.Table));
Xa  = [];

for ii = 1:numel(xa)
    for jj = 1:Shapes.NModal(ii)
        Xa = [Xa,xa(ii)];
    end
end

Shapes.NDof = sum(Shapes.Table);
Shapes.NDim = sum(Shapes.NModal);
Shapes.Ba   = I6(:,Xa);

Shapes.Gamma = [];
Shapes.Kappa = [];

if isempty(varargin)
   DATASET = 1:numel(t);
else 
   DATASET = unique(round(varargin{1}));
   if DATASET(1) <= 0 
      DATASET = numel(t);
   end
end

for ii = DATASET
    
   N = Shapes.Fem.Log.Node{ii};
   R = Shapes.Fem.Log.Rotation{ii};
   S = Shapes.Fem.Log.Stretch{ii};
 
   %P = Shapes.Filter;  
   %Np = P*[N(:,1),N(:,1)*0,N(:,2)];
   
   %[Kf, Gf] = DifferentialGeometry(Shapes,Np);
   
   [Kf, Gf] = ReconstructField(Shapes,N,R,S);

   Shapes.Gamma = [Shapes.Gamma, Gf-1];
   Shapes.Kappa = [Shapes.Kappa, Kf];
end

% SVD decompostion of snapshot reconstructions
[Ur,Sr,~] = svd(full(Shapes.Kappa*Shapes.Kappa.'));
[Uq,Sq,~] = svd(full(Shapes.Gamma*Shapes.Gamma.'));

Er = (diag(Sr).^0.5);
Eq = (diag(Sq).^0.5);

Shapes.PODEnergy{2} = Eq/sum(Eq);
Shapes.PODEnergy{1} = Er/sum(Er);

for ii = 1:10
    PODr = Ur(:,ii);
    Shapes.PODR = [Shapes.PODR,PODr];
end

% gram–schmidt orthogonalization, i.e., int_X y1*y2 ds = 1
Shapes.PODR = gsogpoly(Shapes.PODR,Shapes.Sigma);

for ii = 1:10
    PODq = Uq(:,ii);
    Shapes.PODQ = [Shapes.PODQ, PODq];
end

% gram–schmidt orthogonalization
Shapes.PODQ = gsogpoly(Shapes.PODQ,Shapes.Sigma);

k = 1;
Shapes.POD = [];
for ii = 1:numel(Shapes.NModal)
for jj = 1:Shapes.NModal(ii)
    if ii == 1
        Shapes.POD(:,k) = Shapes.PODR(:,jj);
    else
        Shapes.POD(:,k) = Shapes.PODQ(:,jj);
    end
    k = k+1;
end
end
    
% rebuild shape-function matrix
Shapes.Theta = @(x) ShapeFunction(Shapes,x);

% rebuild shape-function matrix
Shapes = Shapes.rebuild();

end
%--------------------------------------------------------- compute jacobian
function [g, J] = string(Shapes,q)
    
if numel(q) ~= Shapes.NJoint
   error(['Dimension of joint inconstisten with POD matrix. Please ', ...
       'check your input dimensions dim(q).']) 
end    

% ensures robustness for near-zero singularities in some PCC models
q = q(:) + 1e-12;

[g,J] = computeForwardKinematicsFast_mex(q,q*0,... % states
    Shapes.ds,...         % spatial steps
    Shapes.g0(1:3,4),...         % position zero
    Shapes.g0(1:3,1:3),...       % phi zeroclc
    Shapes.Xi0Eval,...    % intrinsic strain vector
    Shapes.ThetaEval,...  % evaluated Theta matrix
    Shapes.Ba);


% gtmp = zeros(4,4,Ns);    % cell(numel(Shapes.Sigma),1);
% Jtmp = zeros(6,Ndim,Ns); % cell(numel(Shapes.Sigma),1);
% 
% for ii = 1:Ns
%     
%    K1 = ForwardKinematicODE(Shapes,s,q,X);
%    K2 = ForwardKinematicODE(Shapes,s+(2/3)*h,q,X + (((2/3)*h)*K1));
%     
%    s = s + h; 
%    X = X + 0.25*h*(K1 + 3*K2);
%    
%    gtmp(:,:,ii) = SE3(X.Phi,X.p);
%    Jtmp(:,:,ii) = Admapinv(X.Phi,X.p)*X.J;
% 
% end
% 
% J = Jtmp;
% g = gtmp;
    
end
%--------------------------------------------------------- compute jacobian
function [dx, Shapes] = flow(Shapes,q,varargin)
if ~isempty(varargin)
    u = varargin{1}(:);
    if numel(varargin)>1
        t = varargin{2};
    end
else
    t = 0;
    u = zeros(Shapes.NInput,1);
end
    
NQ   = Shapes.NJoint;
Q    = q(1:NQ);
dQ   = q(NQ+1:2*NQ);
    
[g,J,eta] = computeForwardKinematicsFast_mex(Q,dQ,... % states
    Shapes.ds,...         % spatial steps
    Shapes.g0(1:3,4),...         % position zero
    Shapes.g0(1:3,1:3),...       % phi zero
    Shapes.Xi0Eval,...    % intrinsic strain vector
    Shapes.ThetaEval,...  % evaluated Theta matrix
    Shapes.Ba);
% 
% [M_,C_,K_,R_,fg_,...
%     gam_,Phi_,J_,dJdt_,Grav,Kin] = computeLagrangianFast_mex(...
%     Q,dQ,...
%     Shapes.ds,...
%     Shapes.g0(1:3,4),...
%     Shapes.g0(1:3,1:3),...
%     Shapes.Xi0Eval,...
%     Shapes.ThetaEval,...
%     Shapes.Ba,...
%     Shapes.Ktt,...
%     Shapes.Mtt,...
%     Shapes.Material.Zeta,...
%     Shapes.Gravity);

[M_,C_,K_,R_,fg_,G_,...
    gam_,Phi_,J_,dJdt_,Grav,Kin] = computeLagrangianFast_mex(...
    Q,dQ,...
    Shapes.ds,...
    Shapes.g0(1:3,4),...
    Shapes.g0(1:3,1:3),...
    Shapes.Xi0Eval,...
    Shapes.ThetaEval,...
    Shapes.MuscleEval,...
    Shapes.Ba,...
    Shapes.Ktt,...
    Shapes.Mtt,...
    Shapes.Material.Zeta,...
    Shapes.Gravity);

% overwrite dynamics
Shapes.Log.EL.M     = M_;
Shapes.Log.EL.C     = C_;
Shapes.Log.EL.R     = R_;
Shapes.Log.EL.K     = K_;
Shapes.Log.EL.fg    = fg_;
Shapes.Log.EL.J     = J_;
Shapes.Log.EL.dJdt  = dJdt_;
Shapes.Log.EL.fc    = 0;
Shapes.Log.EL.dfcdq = 0;

% overwrite energies
Shapes.Log.PH.Kinetic = Kin;
Shapes.Log.PH.Elastic = 0.5*Q.'*K_*Q;
Shapes.Log.PH.Gravity = Grav;
Shapes.Log.PH.dVdq    = K_*Q + fg_;

if isempty(Shapes.Muscle)
   Shapes.Log.EL.G = eye(Shapes.NJoint);
else
   Shapes.Log.EL.G = G_; 
end

% Shapes.Log.EL.G = Shapes.InputMap(Shapes);

% if ~isempty(Shapes.Muscle)
%    P = computeMuscleGroups(Shapes); 
% end

% pre-compute Minverse
Minv = M_\eye(numel(Q));
Shapes.Log.EL.Minv = Minv;

Shapes.Log.t    = t;
Shapes.Log.q    = Q;
Shapes.Log.dq   = dQ;
Shapes.Log.p    = M_*dQ;

Shapes.Log.FK.g    = g;
Shapes.Log.FK.J    = J;
Shapes.Log.FK.eta  = eta;
Shapes.Log.FK.Node = reshape(g(1:3,4,:),3,[]).';
Shapes.Log.FK.L    = sum(sqrt(sum(diff(Shapes.Log.FK.Node).^2,2))) ...
    + Shapes.Length/Shapes.NNode;

Shapes.Log.FK.gam = gam_;
Shapes.Log.FK.Phi = Phi_;

if ~isempty(Shapes.Contact)
   [Shapes.Log.EL.fc,Wc,Wt] = computeContactWrench(Shapes,Q,dQ);

   Nds = Shapes.Log.FK.Node;
   I = ~~sign(sum(abs(Wc),2));
   
   Shapes.Log.Con.Wc = Wc(I,:);    % Wrench normal contact
   Shapes.Log.Con.Wt = Wt(I,:);    % Wrench tangent contact
   Shapes.Log.Con.Node = Nds(I,:);
end

% flow field
dx = [dQ; ...
     Minv*(Shapes.Log.EL.G*u - Shapes.Log.EL.C*dQ - ...
           Shapes.Log.EL.K*Q - Shapes.Log.EL.R*dQ - ...
           Shapes.Log.EL.fg  + Shapes.Log.EL.fc)];
    
end  
%--------------------------------------------------------- compute jacobian
function H = hessian(Shapes,~,varargin)

Minv = Shapes.Log.EL.Minv;
Nq   = Shapes.NJoint;

H                      = zeros(Shapes.NDim,Shapes.NDim);
H(1:Nq,Nq+1:2*Nq)      = eye(Nq);
H(Nq+1:2*Nq,1:Nq)      = -Minv*(Shapes.Log.EL.K + Shapes.Log.EL.dfcdq);
H(Nq+1:2*Nq,Nq+1:2*Nq) = -Minv*(Shapes.Log.EL.R + Shapes.Log.EL.C);

end
%--------------------------------------------------------- compute jacobian
function varargout = FK(Shapes,q,dq)
    
    if nargin < 2
        dq = q*0;
    end
    
    p0           = Shapes.g0(1:3,4);
    [g, J]       = string(Shapes,q);
    varargout{1} = [p0.'; reshape(g(1:3,4,1:end),3,[]).'];
    
    if nargout > 1
        V = zeros(Shapes.NNode,6);
        
        for ii = 1:Shapes.NNode
            Ge = SE3(g(1:3,1:3,ii));
            V(ii,:) = (J(:,:,ii)*dq).';
        end
        
        varargout{2} = V;
    end
    
    if ~isempty(Shapes.gL)
        gg = g(:,:,end)*Shapes.gL;
        p = varargout{1};
        p(end+1,:)   = gg(1:3,4).';
        varargout{1} = p;
    end
    
end
%--------------------------------------------------------- compute jacobian
function xi = strain(Shapes,q)
    
    q  = q(:);
    xi = zeros(Shapes.NNode,6);
    jj = 1;
    for ii = 1:2:Shapes.NNode*2
        xi(jj,:) = (Shapes.Ba*Shapes.ThetaEval(:,:,ii)*q + ...
            Shapes.Xi0Eval(:,:,ii)).';
        jj = jj + 1;
    end
    
end
%------------------------------------------------- estimate Cosserat string
function [q, XI, Nc] = estimateJointSpace(Shapes,Xi,Nds)

Kappa_ = Xi(:,1);
Gamma_ = Xi(:,2);

ki = Shapes.NModal(2);
ei = Shapes.NModal(4);

PODr = Shapes.PODR(:,1:ki);
PODq = Shapes.PODQ(:,1:ei);

if nargin < 3
   W = [1,1];
end

WR = 1;
WQ = 1;

Lam1 = diag(trapz(Shapes.Sigma,PODr.*PODr));
Lam2 = diag(trapz(Shapes.Sigma,PODq.*PODq));

XR = Lam1\trapz(Shapes.Sigma,PODr.*WR.*Kappa_).';
XQ = Lam2\trapz(Shapes.Sigma,PODq.*WQ.*Gamma_).';

q0 = [XR; XQ];

if isempty(Shapes.Filter)
    Shapes = GenerateRadialFilter(Shapes);
end

Nc = Shapes.Filter*Nds;  

k1 = 0.000;
k2 = 0.25;

E = 1;
q = q0;

Itr = 1;
Ids = [round(1:(Shapes.NNode/10):Shapes.NNode),(Shapes.NNode)];

while (E > 1e-3) && (Itr < 150)
    
    [gg, JJ] = Shapes.string(q);
    
    dq = q*0;
    
    for ii = 1:numel(Ids)
        
        g  = gg(:,:,Ids(ii));
        J  = JJ(:,:,Ids(ii));
        gd = SE3(g(1:3,1:3),[Nc(Ids(ii),1);0;Nc(Ids(ii),2)]);
        
        lam1 = Shapes.Kp;
        lam2 = Shapes.Kd;
        
        % conditioner
        Ke = diag([k1,k1,k1,k2,k2,k2])*(ii/numel(Ids)).^4;
        
        Xi = logmapSE3(g\gd);
        Fu = Ke*tmapSE3(Xi)*wedge(Xi);
        
        E = norm(wedge(Xi));
        
        dq = dq + lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
    end
    
    q   = q + dq;
    Itr = Itr + 1;
end

% opts = optimoptions('fminunc','FiniteDifferenceStepSize',1e-6);
% q    = fminunc(@(x) Objective(x,Shapes,Nc(end,:)),q0,opts);
% % 
%     function J = Objective(Q,shp,P)
%         P_ = shp.FK(Q);      
%         J = sum((sum((P - P_(end,[1,3])).^2,2)));  
%         %[~, D] = distance2curve(P,P_(end,[1,3]));
%         %J = sum(D);
%     end

% XR = PODr.'*inv(PODr*PODr.')*Kappa_;
% XQ = PODq.'*inv(PODq*PODq.')*Gamma_;
%q  = q0;
XI = [PODr*XR,PODq*XQ];

end
%------------------------------------------------- estimate Cosserat string
function Xi = recoverStrain(Shapes,fem,varargin)
   
Shapes = Shapes.set('Fem',fem);   
    
if isempty(varargin)
    N = fem.Log.Node{end};
    R = fem.Log.Rotation{end};
    S = fem.Log.Stretch{end};
else
    N = fem.Log.Node{varargin{1}};
    R = fem.Log.Rotation{varargin{1}};
    S = fem.Log.Stretch{varargin{1}};
end

[Curvature, Stretch, ~, ~] = ReconstructField(Shapes,N,R,S);

G = full(Curvature);
U = full(Stretch-1);

Xi = [G,U];
    
end
%----------------------------------------------------- tangent point energy
function [E, R] = tangentPoint(Shapes,q,varargin) 
        
% compute string  tangent
[gc,~] = string(Shapes,q);

if nargin < 3
    PS1 = reshape(gc(1:3,4,1:end),3,[]).';
else
    PS1 = varargin{1};
end

% compute minimal torus
R = zeros(length(gc),1);

for ii = 1:length(gc)

    g = gc(:,:,ii);
    T = g(1:3,1);
    x = g(1:3,4);
    
    %for jj = 1:length(PS1)
        N = (PS1 - x.');
        R2 = diag(N*N.'); 
        R2(abs(R2) <= 1e-6) = max(R2);
        N = N./sqrt(sum(N.^2,2));
        TangX = cross(repmat(T.',length(R2),1),N);
        TangX(isnan(TangX(:,1)),:) = 0;
        T2 = diag(TangX*TangX.');        
        r = T2./R2;
    %end
%     
    R(ii) = sum(r);
end

R = GaussianFilter(R,round(numel(R)*0.01));
E = trapz(Shapes.Sigma,R);

end
%----------------------------------------------------- tangent point energy
end

methods (Access = private)
%---------------------------------------------------------------------- set
function P = ShapeFunction(Shapes,X)

    k  = 1;
    X0 = Shapes.Sigma;
    Pc = cell(Shapes.NDof,1); 
    %X  = zclamp(X,0,Shapes.Length); % make bounded
    X = zclamp(X/Shapes.Length,0,1);

    % construct shape-matrix 
    for jj = 1:6
        for ii = 1:Shapes.NModal(jj)
            
            if jj < 4
                THETA = Shapes.PODR(:,ii);   % angular strains
            else
                THETA = Shapes.PODQ(:,ii);  % linear strains
            end
            % not sure if interp1 is best/fastest option? 
            % maybe inverse lerp?
            Pc{k,1} = interp1(X0,THETA,X);
            k = k + 1;
        end
    end

    P = blkdiag(Pc{:});

end
%---------------------------------------------------------------------- set
function P = IntrinsicFunction(Shapes,X)
    P = Shapes.xia0;
end
%-------------------------------------------------- compute Cosserat string
function Shapes = GenerateRadialFilter(Shapes)
   
PS1 = Shapes.Node0; 
PS2 = Shapes.Fem.get('Center0');
PS3 = Shapes.Fem.get('Node0');
ShapeFnc = Shapes.Fem.get('ShapeFnc');

d = cell(size(PS1,1),1);    

% get closest neighbors with element centers
I = knnsearch(PS1,PS2);

for ii = 1:Shapes.NNode
   el = Shapes.Fem.Element{I(ii)};
   Vsp = PS3(el,:);  % points spanned by element
   Vin = PS1(ii,:);  % point in element
   
   % compute shape function jacobian transformation
   J0 = Vsp.'*ShapeFnc{numel(el)}.dNdxi(:,:,1);
   
   V1 = ((J0)\Vsp.').';     % transform elements to reference config.
   dv = mean(V1,1);
   V1 = V1 - dv;            % pull to origin
   V2 = ((J0)\Vin.').'- dv; % transform linepoint to reference config.

   % recover angle
   th = atan2(V2(2),-V2(1));
   r  = sqrt(dot(V2,V2));
   
   L = [0,0; r*cos(pi+th),-r*sin(pi+th)]; % draw line to point 
   P = [V1;V1(1,:)];   
   [xc, yc] = intersections(P(:,1),P(:,2),L(:,1),L(:,2)); 
   
   if isempty(xc) % outside element
    N = ShapeFnc{numel(el)}.fnc(r*[cos(th),sin(th)]);
   else % inside element    
    % find line intersection on boundary of element
    N = ShapeFnc{numel(el)}.fnc([-xc,yc]);   
   end
   
   % assemble distance filter matrix based on ShpFnc N(s)
   if numel(el) == 3
    d{ii} = [repmat(ii,numel(el),1),[el(2);el(1);el(3)],(N)];
   else
    d{ii} = [repmat(ii,numel(el),1),[el(end-1:-1:1).';el(end)],(N)];   
   end
   
end

d = cell2mat(d); 

P = sparse(d(:,1),d(:,2),d(:,3),Shapes.NNode,Shapes.Fem.NNode);
P = spdiags(1./sum(P,2),0,size(P,1),size(P,1))*P;

Shapes.Filter = P;

end
%-------------------------------------------------- compute Cosserat string
function [Kappa, Gamma, GSE3, Nod] = ReconstructField(Shapes,N,R,S)
    
P   = Shapes.Filter;  
lst = 1:Shapes.Fem.NNode;

Gamma = zeros(Shapes.NNode,1);
Kappa = zeros(Shapes.NNode,1);
Rot = {};
Nod = P*N;
GSE3 = zeros(4,4,Shapes.NNode);

% loop over each node in curve to find geometric strains
for ii = 1:Shapes.NNode
    
    RRe = 0;
    UUe = 0;
    
    W = P(ii,:);
    
    for jj = lst(abs(P(ii,:))>0)
        UUe = UUe + W(jj)*S{jj};
        RRe = RRe + W(jj)*R{jj};
    end
    
    [Ur,~,Vr] = svd(RRe);
    Re = (Ur*Vr.');
     
    Rot{ii,1}      = Re;
    Tangent(ii,:)  = Re(:,1).';
    Gamma(ii,1)    = trace(UUe)/3;
end

for ii = 2:Shapes.NNode-1
    Ni  = Nod(ii,:);
    Nii = Nod(ii+1,:);
    Nip = Nod(ii-1,:);
    
    g  = [Rot{ii,1},    [Ni(1);  Ni(2);  0]; zeros(1,3), 1];
    dg = ([Rot{ii+1,1}, [Nii(1); Nii(2); 0]; zeros(1,3), 1] ...
        - [Rot{ii-1,1}, [Nip(1); Nip(2); 0]; zeros(1,3), 1])/(2*Shapes.ds);
    
    GSE3(:,:,ii) = g;
    XI = g\dg;
    
    Kappa(ii,1) = XI(1,2);
    Gamma(ii,1) = XI(1,4);
end


% for ii = 1:Shapes.NNode
%     if ii == 1
%         Kappa(ii) = 0;
%     elseif ii == Shapes.NNode
%         Kappa(ii) = Kappa(ii-1);
%     else
%         t1 = (Tangent(ii,:) + Tangent(ii-1,:)); t1 = t1/norm(t1);
%         t2 = (Tangent(ii,:) + Tangent(ii+1,:)); t2 = t2/norm(t2);
%         dir = sign(dot(so3(t2)*t1(:),[0,0,1])); 
%         angle = real(2*dir*acos(dot(t2,t1)));
%         
%         % differential geometric on discretized curve 
%         % recover the approximate curvature
%         Kappa(ii,1) = angle/(norm(Nod(ii+1,:) - Nod(ii,:)) + ...
%             norm(Nod(ii,:) - Nod(ii-1,:)));
%     end   
% 
% end


Kappa(1,1) = Kappa(2,1);
Gamma(1,1) = Gamma(2,1);
Kappa(end,1) = Kappa(end-1,1);
Gamma(end,1) = Gamma(end-1,1);

if ~isempty(Shapes.FilterRadius)
    Rf = Shapes.FilterRadius;
    if numel(Rf) == 1
        Kappa = GaussianFilter(Kappa,round(Rf));
        Gamma = GaussianFilter(Gamma,round(Rf));
    else 
        Kappa = smoothdata(Kappa,'gaussian',round(Rf(1)));%sgolayfilt(Kappa,3,round(Rf(1)));
        Gamma = smoothdata(Gamma,'gaussian',round(Rf(2)));%sgolayfilt(Gamma,3,round(Rf(2)));
    end
end

% subplot(2,1,1);
% plot(Gamma); hold on
% 
% subplot(2,1,2);
% plot(Kappa); hold on

end
%------------------------------------- compute differential geometry curve
function [Kappa, Gamma] = DifferentialGeometry(Shapes,Node)
% http://page.math.tu-berlin.de/~bobenko/Lehre/Skripte/DDG_Lectures.pdf
N = Shapes.NNode;
T = zeros(N,3); 
%dgam = zeros(N,3); 
%phi  = zeros(N,1);

[~, Fy] = gradient(Node);   
dgam = [Fy(:,1),Fy(:,2),Fy(:,3)]; 

[~, Fy] = gradient(Shapes.Node0);   
dgam0 = [Fy(:,1),Fy(:,2)]; 

dl  = sqrt(sum(dgam.^2,2));
dl0 = sqrt(sum(dgam0.^2,2));

Gamma = dl./dl0;

% compute tangents
for ii = 1:N
    if ii < Shapes.NNode
        dgam(ii,:) = Node(ii + 1,:) -  Node(ii,:);
    else
        dgam(ii,:) = Node(ii,:) - Node(ii - 1,:);
    end
    
    T(ii,:) = dgam(ii,:)/norm(dgam(ii,:));
    
end

I = null(round(T(1,:)));
Normal = I(:,1);
    
% compute curvature
for ii = 2:Shapes.NNode-1
    
    %t1 = T(ii - 1,:);
    %t2 = T(ii,:);
%     
    t1 = (T(ii,:) + T(ii-1,:)); t1 = t1/norm(t1);
    t2 = (T(ii,:) + T(ii+1,:)); t2 = t2/norm(t2);
    
    dir = sign(dot(so3(t2)*t1(:),Normal));
    angle = real(2*dir*acos(dot(t2,t1)));
    
    %dir = sign(dot(so3(t2)*t1(:),Normal));
    %phi = real(2*dir*acos(dot(t2,t1)));
    
    Kappa(ii,1) = angle/(norm(Node(ii+1,:) - Node(ii,:)) + ...
        norm(Node(ii,:) - Node(ii-1,:)));
end

Kappa(1,:) = Kappa(2,:) + Shapes.ds*(Kappa(3,:) - Kappa(2,:));
Kappa(N,:) = Kappa(end-1,:) + Shapes.ds*(Kappa(end,:) - Kappa(end-1,:));

%Kappa    = GaussianFilter(Kappa,round(Shapes.NNode/1000));
%Gamma    = GaussianFilter(Gamma,round(Shapes.NNode/1000));

end
%--------------------------------------- forwards integration of kinematics
function dg = ForwardODE(Shapes,s,g,q)
    
%compute strain field    
xi = Shapes.Ba*Shapes.Theta(s)*q + Shapes.xia0(:);

Kap = xi(1:3);  % get curvature-torsion
Gam = xi(4:6);  % get stretch-shear
Q   = g(1:4);   % get quaternions

R = Quat2Rot(Q); % thanks for your inspiring work Frederic Boyer ;)
A = StrainMap(R*Kap(:));

dg      = zeros(7,1);
dg(1:4) = ((2*norm(Q)))\A*Q;
dg(5:7) = R*Gam(:);

end    
%--------------------------------------- forwards integration of kinematics
function Y = ForwardKinematicODE(Shapes,s,q,X)

% recover variables from struct
Y    = X;
Phi_ = X.Phi;
p_   = X.p;
    
%compute strain field    
xi = Shapes.Ba*Shapes.Theta(s)*q + Shapes.xia0(:);

% construct geometric vectors
Kap = xi(1:3);  % get curvature-torsion
Gam = xi(4:6);  % get stretch-shear

% build forward kin - position
Y.p   = Phi_*Gam;
Y.Phi = Phi_*isomSO3(Kap);
Y.J   = Admap(Phi_,p_)*Shapes.Ba*Shapes.Theta(s);

end
%--------------------------------------- forwards integration of kinematics
function Shapes = BuildInertia(Shapes)
    
N  = round(Shapes.NNode/1.5);
x0 = linspace(Shapes.Sdf.BdBox(1),Shapes.Sdf.BdBox(2),N);
y0 = linspace(Shapes.Sdf.BdBox(3),Shapes.Sdf.BdBox(4),N);
[X0,Y0] = meshgrid(x0,y0);

% get tangent-sub volume
dv = (x0(2) - x0(1))*(y0(2) - y0(1));

% generate image from cross-section
D   = Shapes.Sdf.eval([X0(:),Y0(:)]);
rho = (D(:,end)<1e-5);

I0 = reshape(rho,[N,N]);

% https://ocw.mit.edu/courses/aeronautics-and-astronautics/
% 16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec26.pdf
x0 = x0 - Shapes.Center(1);
y0 = y0 - Shapes.Center(2);
X0 = X0 - Shapes.Center(1);
Y0 = Y0 - Shapes.Center(2);

% evaluate slice volume
Shapes.Att = sum(sum(I0*dv));

% evaluate 2nd-moment inertia
Jxx = trapz(y0,trapz(x0,(Y0.^2).*I0,2))/Shapes.Att;
Jyy = trapz(y0,trapz(x0,(X0.^2).*I0,2))/Shapes.Att;
Jzz = trapz(y0,trapz(x0,(X0.^2 + Y0.^2).*I0,2))/Shapes.Att;

Jxy = trapz(y0,trapz(x0,(X0.*Y0).*I0,2))/Shapes.Att;
Jxz = 0*trapz(y0,trapz(x0,(X0.*Y0).*I0,2))/Shapes.Att;
Jyz = 0*trapz(y0,trapz(x0,(X0.*Y0).*I0,2))/Shapes.Att;

P = eye(3);
P = [P(:,3),P(:,1),P(:,2)];

Shapes.Jtt = P.'*[Jxx,Jxy,Jxz;Jxy,Jyy,Jyz;Jxz,Jyz,Jzz]*P;
Shapes.Mtt = Shapes.Material.Rho*blkdiag(Shapes.Jtt,Shapes.Att*eye(3));

end
%--------------------------------------- forwards integration of kinematics
function [Fr,Wc,Wt]= computeContactWrench(Shapes,Q0,dQ0)
 
Wc = zeros(3,Shapes.NNode);
Wt = zeros(3,Shapes.NNode);    
Fr = zeros(Shapes.NJoint,1);

if ~isempty(Shapes.Log.FK)
    [g, J] = Shapes.string(Q0);
else
    g   = Shapes.Log.FK.g;
    J   = Shapes.Log.FK.J;
end

Nqu  = 6; 
Emod = Shapes.Material.getModulus();
Rmod = Shapes.Material.getContactReaction();
Cmod = Shapes.Material.getContactFriction();
Eps  = Shapes.ContactDistance;

% generate vector of radii
if Shapes.VolumetricContact
    R = ones(Shapes.NNode,1)*Shapes.TubeRadiusA;
    X = 1:Shapes.NNode;
    
    if ~isempty(Shapes.TubeRamp)
        y = Shapes.TubeRamp;
        if numel(y) == 1
            R = R(1,1)*(1-y*(X/Shapes.NNode));
        else
            S = linspace(1,Shapes.NNode,numel(y));
            R = interp1(S,R(1,1)*y,X);
        end
    end
end

Pos = reshape(g(1:3,4,1:end),3,[]).';
Id = Shapes.Contact.intersect(Pos,Shapes.TubeRadiusA*2); 
Set = 1:Shapes.NNode;
DetectionSet = Set(Id);

for ii = DetectionSet

    % generate sample set
    if Shapes.VolumetricContact
        XY = SampleRing(R(ii),Nqu,g(:,:,ii));
    else
        XY = g(1:3,4,ii).';
    end
    
    D = Shapes.Contact.eval(XY); 
    D = D(:,end);
    
    ConU  = D(D<Eps);
    ConXY = XY(D<Eps,:);
    JJ = J(4:6,:,ii);
    
    if ~isempty(ConXY)
        [T_,~,B_] = Shapes.Contact.normal(g(1:3,4,ii).');
        [~,N]  = Shapes.Contact.normal(ConXY);

        Vt = JJ*dQ0;
        VT = Vt./(sqrt(sum((Vt.^2))) + 1e-3);
        FT = (-T_.*dot(VT.',T_.').').';
        FB = (-B_.*dot(VT.',B_.').').';
    end
    
    Rw = zeros(3,1);
    Rf = zeros(3,1);
    
    for jj = 1:sum(D<Eps)     
        U = (N(jj,:).')*(ConU(jj));
        Rw = Rw - (1/Nqu)*Emod*Rmod*U(:);
    end

    if ~isempty(ConXY)
        Rf = Cmod*(dot(Rw,FT)).*FT + Cmod*(dot(Rw,FB)).*FB ;
    end

    % project wrench onto conf. space
    Fr = JJ.'*(Rw + Rf);
    
    Wc(:,ii) = Rw;
    Wt(:,ii) = Rf;
end

Shapes.Log.EL.fc = Fr;
Wc = Wc.';
Wt = Wt.';

% -------- 
    function XY = SampleRing(R,Quality,Gd)
        th = linspace(0,2*pi,Quality+1);
        th = th(1:end-1);
        N0 = th(:)*0;
        Circ0 = [N0,R*cos(th(:)),R*sin(th(:)),N0 + 1];
        Circ = (Gd*Circ0.').';
        XY = Circ(:,1:3);
    end
% --------
end
%--------------------------------------------- isomorphism from R3 to so(3)
function Kc = computeContactJacobian(Shapes,Q0,dQ0)
   
n = Shapes.NJoint;    
epsilon = 1e-3;
delta   = epsilon*eye(n);

[Fc0] = computeContactWrench(Shapes,Q0,dQ0);
Kc  = zeros(n);
Dc  = zeros(n);

% finite difference for tau(q(t),.)
for ii = 1:n
    Q  = Q0 + delta(:,ii);
    [Fc_] = computeContactWrench(Shapes,Q,dQ0);
    Kc(:,ii) = (Fc_ - Fc0)/epsilon;
end

end
%--------------------------------------------- isomorphism from R3 to so(3)
function [P,F] = computeMuscleGroups(Shapes,s)
 
M = numel(Shapes.Muscle);
P = zeros(numel(s),3,M);
F = zeros(numel(s),3,M);
r = Shapes.TubeRadiusA;
if isempty(Shapes.TubeRamp)
    R = @(x) r;
else
    R = @(x) r * (1 -  Shapes.TubeRamp * x);
end

Xi03 = Shapes.xia0(4:6);

for ii = 1:M
    fnc = Shapes.Muscle{ii};
    P0  = R(s).*fnc(s);
    P(:,:,ii) = P0.' + (s.*(Xi03(:))*Shapes.Length).';
    F(:,:,ii) = P0.';
end

end

end
end

%--------------------------------------------- isomorphism from R3 to so(3)
function y = isomSO3(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; x3, 0, -x1; -x2, x1, 0];
end
%--------------------------------------------- isomorphism from R3 to so(3)
function y = voightextraction(X)
y(1,1) = X(2,3);
y(2,1) = 2*X(4,4);
y(3,1) = 2*X(4,4);
y(4,1) = X(1,1);
y(5,1) = X(2,1);
y(6,1) = X(3,1);
end
% %------------------------------------------------------- optimal sphere fit
% function r = torusSolve(gp,x)
% % x  := [Nx3] of points
% % gp := [4x4] SE(3) matrix of tangent-point orgin 
% R0 = fliplr(gp(1:3,1:3));
% x0 = gp(1:3,4);
% %x = 
% x = x + rand(length(x),3)*1e-5;
% dx = x.' - x0;
% 
% X1 = (R0(1,:)*(dx)).^2;
% X2 = (R0(2,:)*(dx)).^2;
% X3 = (R0(3,:)*(dx)).^2;
% 
% r = sqrt((((X1 + X2 + X3).^2) ...
%      ./(4*(X1 + X2))));
%  
% r(isnan(r)) = 0;
%  
% end