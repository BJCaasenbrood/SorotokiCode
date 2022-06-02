classdef Shapes

    properties (Access = public)
        Fem;   
        Sdf;
        Material;
        Sigma;
        L0;
        Ba;
        ds;
        xia0;
        Node;
        Node0;
        
        NModal;
        NDim;
        NNode;
        
        POD;
        PODR;
        PODQ;
        Theta;
        Xi0;
        
        Rho;
        Zeta;
        E;
        Nu;
        
        Mtt;
        Dtt;
        Ktt; Ktt0;
        HypA; HypB;
        Jtt;
        Att;
        Gvec;
    end
    
    properties (Access = private)
        Table;
        NDof;
        Center;
        
        Rotation;
        Gamma;
        Kappa;
        
        PODEnergy;
        
        Phi0;
        g0;
        Filter;
        Quality;
        
        ThetaEval;
        Xi0Eval;
    end
   
%--------------------------------------------------------------------------
methods  
%----------------------------------------------- MODAL SHAPE RECONSTRUCTION
function obj = Shapes(Input,NModal,varargin) 
    
    obj.NModal = NModal;
    obj.Table  = double(NModal > 0); 
    obj.NDof   = sum(obj.Table);
    obj.xia0   = [0,0,0,1,0,0].';
    
    if isa(Input,'Fem')
        obj.Fem    = Input;
        obj.NNode  = 30;        
        obj.Rho    = Input.Material.Rho;
        obj.Zeta   = Input.Material.Zeta;
        obj.Gvec   = [0;0;9.81e3];
        
    elseif isa(Input,'double')
        obj.NNode = size(Input,1);
        obj.PODQ = Input;
        obj.PODR = Input;
        
        obj.L0    = 1;
        obj.Sigma = linspace(0,1,obj.NNode); 
        obj.ds    = obj.L0/(obj.NNode);

        obj.PODEnergy{2} = ones(size(Input,2),1);
        obj.PODEnergy{1} = ones(size(Input,2),1);
        
        obj.Rho  = 1090e-12;
        obj.Zeta = .15;
        obj.Gvec = [0;0;9.81e3];
    end

    % cross-section SDF
    obj.Sdf    = sCircle(5);
    obj.Center = zeros(3,1);
    obj.Material = NeoHookeanMaterial(1,0.33);

    obj.E  = 5;
    obj.Nu = 0.33;
    obj.g0 = SE3(eye(3),zeros(3,1));
    obj.HypA = 0;
    obj.HypB = 0;
    
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
    
    %obj.g0 = [1,0,0,0,0,0,0];
    %obj.g0 = SE3(eye(3),zeros(3,1));
    
%     for ii = 1:2:length(varargin)
%         obj.(varargin{ii}) = varargin{ii+1};
%     end
%     
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
%---------------------------------------------------------------- show mesh
function Shapes = show(Shapes,varargin)
    
if nargin<2, Request = -1; 
else, Request = varargin{1}; 
end

figure(101);

switch(Request)
    case('Base'),  Request = 'POD';
    case('POD'),   Request = 'POD';
    otherwise,     Request = 'Curve';
end

if strcmp(Request,'POD') 
    
%     for jj = 1:2
%         
%         % pick shape-matrix
%         if jj == 1, theta = Shapes.PODR; 
%         else, theta = Shapes.PODQ;
%         end
%             
%         k   = 1;
%         leg = cell(1);
%         % loop through selected modes
%         for ii = 1:Shapes.NModal(jj)
%             subplot(Shapes.NDof,2,2*jj - 1)
%             plot(Shapes.Sigma/Shapes.L0,theta(:,ii),...
%                 'Color',col(k),'linewidth',3); 
%             hold on; grid on;
%             k = k +1;
%             leg{ii} = ['\theta_',num2str(ii)];
%         end
%         
%         legend(leg{:},'Fontsize',14,'Orientation','Horizontal');        
%         
%         A  = axis; Ay = 1.2*max(abs(A(3:4)));
%         
%         subplot(Shapes.NDof,2,2*jj);
%         En = (Shapes.PODEnergy{jj})/sum(Shapes.PODEnergy{jj});
%         plot(En,'-o','LineW',3); hold on;
%         Enp = cumsum(En);
%         plot(Enp,'--','LineW',3); 
% 
%     end

end

end
%------------------------------------------------------------ set reference
function Shapes = reference(Shapes,varargin)
X0 = varargin{1};
XL = varargin{2};

% build Cosserat curve in reference
Shapes.Node0 = [linspace(X0(1),XL(1),Shapes.NNode).', ...
                linspace(X0(2),XL(2),Shapes.NNode).'];
             
%B = boxhull(Shapes.Node0);                 
Shapes.L0 = sqrt((XL(1)-X0(1))^2 + (XL(2)-X0(2))^2);  

% build discretization of curve domain
Shapes.Sigma = linspace(0,Shapes.L0,Shapes.NNode);    
Shapes.ds    = Shapes.L0/(Shapes.NNode);

% finds associated nodes from Fem mesh.
if isempty(Shapes.Fem)
    Shapes = GenerateRadialFilter(Shapes);
end

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

Shapes.NDim  = sum(Shapes.NModal);
Shapes.Ba    = I6(:,Xa);
Shapes.Sigma = linspace(0,1,Shapes.NNode);
Shapes.ds    = Shapes.L0/(Shapes.NNode);

Shapes = BuildInertia(Shapes);

JJ     = Shapes.Mtt/Shapes.Material.Rho;
[SS,KK] = Shapes.Material.PiollaStress(eye(3));
KK = 2*diag(voightextraction(KK));
%KK2 = LinearStiffnessTensor(Shapes);
% KQR = KK(4:6,4:6); 
% KQG = KK(1:3,1:3);
% 
% KK = blkdiag(KQR,KQG);
Shapes.Ktt  = JJ*KK;%diag(diag(JJ))*diag(diag(KK));
%Shapes.Mtt  = diag(diag(Shapes.Mtt));
%Shapes.Ktt0 = Shapes.Ktt;
Shapes.Dtt  = Shapes.Material.Zeta*Shapes.Ktt;
  
if ~isempty(Shapes.Node0)
    Shapes = GenerateRadialFilter(Shapes);
end

if ~isempty(Shapes.PODR) || ~isempty(Shapes.PODQ)
    
    if length(Shapes.PODQ) ~= Shapes.NNode
        X = linspace(0,Shapes.L0,length(Shapes.PODQ));
        Shapes.PODR = interp1(X,Shapes.PODR,Shapes.Sigma*Shapes.L0);
        Shapes.PODQ = interp1(X,Shapes.PODQ,Shapes.Sigma*Shapes.L0);
    end
    
%     % ensure orthonormality
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
s   = sort([Shapes.Sigma*Shapes.L0,...
            Shapes.Sigma*Shapes.L0+(2/3)*Shapes.ds]);

[nx,ny]     = size(FncT(0));
Shapes.ThetaEval = zeros(nx,ny,numel(s));
Shapes.Xi0Eval   = zeros(6,1,numel(s));

for ii = 1:numel(s)
    Shapes.ThetaEval(:,:,ii) = FncT(s(ii));
    Shapes.Xi0Eval(:,1,ii)   = FncX(s(ii));
end
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
    
if numel(q) ~= Shapes.NDim
   error(['Dimension of joint inconstisten with POD matrix. Please ', ...
       'check your input dimensions dim(q).']) 
end    

% ensures robustness for near-zero singularities in some PCC models
q = q(:) + 1e-6;

[g,J] = computeForwardKinematicsFast_mex(q,... % states
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
function varargout = FK(Shapes,q,dq)
    
    if nargin < 2
        dq = q*0;
    end
    
    p0 = Shapes.g0(1:3,4);
    [g, J] = string(Shapes,q);
    varargout{1} = [p0.';reshape(g(1:3,4,1:end-1),3,[]).'];
    
    if nargout > 1
        V = zeros(Shapes.NNode,6);
        for ii = 1:Shapes.NNode
            V(ii,:) = (J(:,:,ii)*dq).';
        end
        varargout{2} = V;
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
function [q,XI] = estimateJointSpace(Shapes,Xi,W)

Kappa_ = Xi(:,1);
Gamma_ = Xi(:,2);

ki = Shapes.NModal(2);
ei = Shapes.NModal(4);

PODr = Shapes.PODR(:,1:ki);
PODq = Shapes.PODQ(:,1:ei);

if nargin < 3
   W = [1,1];
end

WR = (W(:,1));
WQ = (W(:,2) );
% 
Lam1 = diag(trapz(Shapes.Sigma,PODr.*PODr));
Lam2 = diag(trapz(Shapes.Sigma,PODq.*PODq));
% % 
XR = Lam1\trapz(Shapes.Sigma,PODr.*WR.*Kappa_).';
XQ = Lam2\trapz(Shapes.Sigma,PODq.*WQ.*Gamma_).';
% 
% XI = [PODr*XR,PODq*XQ];

% figure(106); clf;
% subplot(2,1,1);
% plot(Kappa_); hold on;
% plot(PODr*XR);
% 
% subplot(2,1,2);


% XR = (diag(WR)*PODr).'*inv(PODr*diag(WR)*(PODr.'))*Kappa_;
% XQ = (diag(WQ)*PODq).'*inv(PODq*diag(WQ)*(PODq.'))*Gamma_;

% plot(Gamma_); hold on;
% plot(PODq*XQ);

% q0 = zeros(sum(Shapes.NModal),1);
% 
% q = fminunc(@(x) Objective(x,Shapes,N_),q0);
% 
%     function J = Objective(Q,shp,P)
%        
%         P_ = shp.FK(Q);
%         
%         J = sum((sum((P - P_(2:end,:)).^2,2)));
%         
%     end

% XR = PODr.'*inv(PODr*PODr.')*Kappa_;
% XQ = PODq.'*inv(PODq*PODq.')*Gamma_;
q  = [XR; XQ];
XI = [PODr*XR,PODq*XQ];

end
%------------------------------------------------- estimate Cosserat string
function Xi = recoverStrain(Shapes,fem,varargin)
   
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
    %X  = zclamp(X,0,Shapes.L0); % make bounded
    X = zclamp(X/Shapes.L0,0,1);

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
%    
%     k  = 1;
%     X0 = Shapes.Sigma;
%     X = zclamp(X/Shapes.L0,0,1);  
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
   d{ii} = [repmat(ii,numel(el),1),[el(2);el(1);el(3)],(N)];
   
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
     
    Rot{ii,1} = Re;
    Tangent(ii,:)  = Re(:,1).';
    Gamma(ii,1)    = trace(UUe)/3;
end

for ii = 2:Shapes.NNode-1
    Ni  = Nod(ii,:);
    Nii = Nod(ii+1,:);
    Nip = Nod(ii-1,:);
    
    g  = [Rot{ii,1},    [Ni(1); Ni(2);0];   zeros(1,3), 1];
    dg = ([Rot{ii+1,1}, [Nii(1); Nii(2);0]; zeros(1,3), 1] ...
        - [Rot{ii-1,1}, [Nip(1); Nip(2);0]; zeros(1,3), 1])/(2*Shapes.ds);
    
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

Kappa = GaussianFilter(Kappa,round(Shapes.NNode/150));
Gamma = GaussianFilter(Gamma,round(Shapes.NNode/150));

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
function [dp,dPhi,dJ] = ForwardKinematicODE2(Shapes,s,q,Phi_,p_)
    
%compute strain field    
xi = Shapes.Ba*Shapes.Theta(s)*q + Shapes.xia0(:);

% construct geometric vectors
Kap = xi(1:3);  % get curvature-torsion
Gam = xi(4:6);  % get stretch-shear

% build forward kin - position
dp   = Phi_*Gam;
dPhi = Phi_*isomSO3(Kap);
dJ   = Admap(Phi_,p_)*Shapes.Ba*Shapes.Theta(s)*q;

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
Shapes.Mtt = Shapes.Rho*blkdiag(Shapes.Jtt,Shapes.Att*eye(3));

end
%--------------------------------------- forwards integration of kinematics
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
%----------------------------------------------------------- strain mapping
function Ktt = LinearStiffnessTensor(Shapes)
E0 = Shapes.Material.E;
G0 = (E0)/(2*(1+Shapes.Material.Nu));

QE  = diag([E0,G0,G0]);
QR  = diag([G0,E0,E0]);
Ktt = blkdiag(QR,QE);
end
%------------------------------------------------------- optimal sphere fit
function r = torusSolve(gp,x)
% x  := [Nx3] of points
% gp := [4x4] SE(3) matrix of tangent-point orgin 
R0 = fliplr(gp(1:3,1:3));
x0 = gp(1:3,4);
%x = 
x = x + rand(length(x),3)*1e-5;
dx = x.' - x0;

X1 = (R0(1,:)*(dx)).^2;
X2 = (R0(2,:)*(dx)).^2;
X3 = (R0(3,:)*(dx)).^2;

r = sqrt((((X1 + X2 + X3).^2) ...
     ./(4*(X1 + X2))));
 
r(isnan(r)) = 0;
 
end