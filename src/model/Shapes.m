classdef Shapes

    properties (Access = public)
        Fem;   
        Sdf;
        Sigma;
        L0;
        Ba;
        ds;
        xia0;
        
        NDof;
        NModal;
        NDim;
        NNode;
        
        POD;
        Theta;
        
        Rho;
        Zeta;
        E;
        Nu;
        
        Mtt;
        Dtt;
        Ktt; Ktt0;
        Jtt;
        Att;
        
    end
    
    properties (Access = private)
        Table;
        Node;
        Node0;
        Center;
        
        Rotation;
        Gamma;
        Kappa;
        
        PODR;
        PODQ;
        PODEnergy;
        
        Phi0;
        g0;
        Filter;
        Quality;
    end
   
%--------------------------------------------------------------------------
methods  
%----------------------------------------------- MODAL SHAPE RECONSTRUCTION
function obj = Shapes(Fem,Table,varargin) 
    
    obj.Table  = Table; 
    obj.Fem    = Fem;  
    obj.NNode  = 30;
    obj.NModal = 2;
    obj.NDof   = sum(Table);
    
    obj.xia0 = [0,0,0,1,0,0].';
    
    % cross-section SDF
    obj.Sdf = sCircle(5);
    obj.Center = zeros(3,1);
    obj.Rho  = Fem.Material.Rho;
    obj.Zeta = Fem.Material.Zeta;
    obj.E    = 5;
    obj.Nu   = 0.33;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    if ~isempty(obj.Fem.get('Output'))
        out = obj.Fem.get('Output');
        Nd0 = obj.Fem.get('Node0');
        
        BdBox = boxhull(Nd0(out(:,1),:));
        obj   = reference(obj,[BdBox(1),BdBox(3)],...
                              [BdBox(2),BdBox(2)]);
    end
    
    obj.g0 = [1,0,0,0,BdBox(1),0,BdBox(3)];
    
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
    
if nargin<2
    Request = -1; 
else
    Request = varargin{1}; 
end

figure(101);

switch(Request)
    case('Base'),  Request = 'POD';
    case('POD'),   Request = 'POD';
    otherwise,     Request = 'POD';
end
    
if strcmp(Request,'POD') 
    
    for jj = 1:2
        
        % pick shape-matrix
        if jj == 1, theta = Shapes.POD; 
        else, theta = Shapes.PODQ;
        end
            
        k = 1;
        leg = cell(1);
        % loop through selected modes
        for ii = 1:Shapes.NModal(jj)
            subplot(Shapes.NDof,2,2*jj - 1)
            plot(Shapes.Sigma/Shapes.L0,theta(:,ii),...
                'Color',col(k),'linewidth',3); 
            hold on; grid on;
            k = k +1;
            leg{ii} = ['\theta_',num2str(ii)];
        end
        
        legend(leg{:},'Fontsize',14,'Orientation','Horizontal');        
        
        A  = axis; Ay = 1.2*max(abs(A(3:4)));
        axis([0 1 -Ay Ay]);
        
        subplot(Shapes.NDof,2,2*jj);
        En = (Shapes.PODEnergy{jj})/sum(Shapes.PODEnergy{jj});
        plot(En,'-o','LineW',3); hold on;
        Enp = cumsum(En);
        plot(Enp,'--','LineW',3); 
        axis([1 10 0 max(En)*1.2]);
        
        
    end

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
Shapes.ds = Shapes.L0/(Shapes.NNode);

% finds associated nodes from Fem mesh.
Shapes = GenerateRadialFilter(Shapes);

end   
%------------------------------------------------------------ set reference
function Shapes = rebuild(Shapes,varargin)
    
for ii = 1:2:length(varargin)
    Shapes.(varargin{ii}) = varargin{ii+1};
end

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

Shapes = BuildInertia(Shapes);
Shapes.Dtt = Shapes.Zeta*Shapes.Mtt;
Shapes.Ktt = LinearStiffnessTensor(Shapes);
Shapes.Ktt0 = Shapes.Ktt;
    
Shapes = GenerateRadialFilter(Shapes);

end 
%--------------------------------------------------------------------------
function Shapes = reconstruct(Shapes)
    
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

for ii = 1:numel(t)
    
   N = Shapes.Fem.Log.Node{ii};
   R = Shapes.Fem.Log.Rotation{ii};
   S = Shapes.Fem.Log.Stretch{ii};
 
   %P = Shapes.Filter;  
   %Np = P*[N(:,1),N(:,1)*0,N(:,2)];
   
   %[Kf, Gf] = DifferentialGeometry(Shapes,Np);
   
   [~, ~, Kf, Gf] = ReconstructField(Shapes,N,R,S);

   Shapes.Gamma = [Shapes.Gamma - 1, Gf];
   Shapes.Kappa = [Shapes.Kappa, Kf];
end

% SVD decompostion of snapshot reconstructions
[Ur,Sr,~] = svd(full(Shapes.Kappa*Shapes.Kappa.'));%/numel(t));
[Uq,Sq,~] = svd(full(Shapes.Gamma*Shapes.Gamma.'));%/numel(t));

Er = (diag(Sr).^0.5);
Eq = (diag(Sq).^0.5);

Shapes.PODEnergy{2} = Eq;
Shapes.PODEnergy{1} = Er;

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

end
%------------------------------------------------------------------ fitting
function [Shapes, P] = fitPOD(Shapes)
    
    x0 = zeros(Shapes.NDim,1);    
    X  = cell(length(Shapes.posData)-1,1);
        
    % setting optimization settings
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
        'FiniteDifferenceStepSize',Shapes.DiffStepSize);
    
    for jj = 1:Shapes.NDof
        Shapes.xiData{jj} = zeros(Shapes.Quality,1);
    end
    
    for ii = 2:4:length(Shapes.posData)
        
        fun = @(x) Objective(x,Shapes,ii);

        % solve optimization problem
        x = fmincon(fun,x0,[],[],[],[],[],[],[],options);
        X{ii-1} = x;
        
        % overwrite initial conditions
        x0 = x;
        
        cla; 
        [P,~,S] = string(Shapes,x);
        plot(P(:,1),P(:,3),'-'); axis equal; hold on;
        N = Shapes.posData{ii};
        plot(N(:,1),N(:,2),'.'); axis equal; hold on;
        pause(0.5);
        
        ee = zeros(Shapes.NDof,length(S));
        for jj = 1:length(S)
            ee(:,jj) = Shapes.Phi(S(jj))*x + Shapes.Ba.'*Shapes.xia0(:);
        end
        
        for jj = 1:Shapes.NDof
            Shapes.xiData{jj} = vappend(Shapes.xiData{jj},ee(jj,:).',2);
        end
    end
    
    XI = Shapes.xiData{1};
    C = XI*XI.';
    
    [Shapes.POD,E,~] = svd(C);
    Shapes.PODEnergy = diag(E);

    function f = Objective(x,Shapes,id)
       p = string(Shapes,x);
       
       if size(Shapes.posData{id},2) == 2
           p = [p(:,1),p(:,3)];
       end
       
       [~,dist] = distance2curve(p,Shapes.posData{id},'linear'); 
       f1 = sum(abs(dist));
       
       Gl = Shapes.posData{id}(end,:);
       Pl = p(end,:);
       f2 = norm(abs(Gl - Pl));
       
       f = f1 + 0.01*f2;
    end
    
end
%-------------------------------------------------- compute Cosserat string
function [p, R] = string(Shapes,q)

if numel(q) ~= Shapes.NDim
   error(['Dimension of joint inconstisten with POD matrix. Please ', ...
       'check your input dimensions dim(q).']) 
end
    
% forward (explicit) integration of kinematics    
[~,g] = ode23t(@(t,x) ForwardODE(Shapes,t,x,q(:)),...
    Shapes.Sigma,Shapes.g0);

q = g(:,1:4);               % quaternions
p = [g(:,5),g(:,6),g(:,7)]; % position vector

R = {};
for ii = 1:Shapes.NNode
    R{ii,1} = quat2rot(q(ii,:));
end
    
end
%------------------------------------------------- estimate Cosserat string
function q = estimateJointSpace(Shapes,fem)

P = Shapes.Filter;  
N = fem.Log.Node{end};
R = fem.Log.Rotation{end};
S = fem.Log.Stretch{end};

[~, ~, K, G] = ReconstructField(Shapes,N,R,S);
%
% 
%N_    = P*[N(:,1),N(:,1)*0,N(:,2)];
%[K,~] = DifferentialGeometry(Shapes,N_);
% K = 0.5*(K1 + K2);
  
Kappa_ = full(-K);
Gamma_ = full(G-1);

ki = Shapes.NModal(1);
ei = Shapes.NModal(2);

PODr = Shapes.POD(:,1:ki);
PODq = Shapes.POD(:,ki+1:ki+ei);

%XR = trapz(Shapes.Sigma,PODr.*Kappa_);
%XQ = trapz(Shapes.Sigma,PODq.*(Gamma_));
%q  = [XR(:); XQ(:)];

% figure(106);
% subplot(2,1,1);
% plot(Kappa_); hold on;
% plot(PODr*XR.');
% 
% subplot(2,1,2);
% plot(Gamma_-1); hold on;
% plot(PODq*XQ.');
% 

% q = fminunc(@(x) Objective(x,Shapes,N_),q0);
% 
%     function J = Objective(Q,shp,P)
%        
%         P_ = shp.string(Q);
%         
%         J = sum((sum((P - P_).^2,2)));
%         
%     end

XR = PODr.'*inv(PODr*PODr.' + 1e-4*eye(Shapes.NNode))*Kappa_;
XQ = PODq.'*inv(PODq*PODq.' + 1e-4*eye(Shapes.NNode))*Gamma_;

q = [XR;XQ];

end

end

methods (Access = private)
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
   Vin = PS1(ii,:); % point in element
   
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
%---------------------------------------------------------------------- set
function P = ShapeFunction(Shapes,X)

    Pc  = cell(Shapes.NDof,1); 
    
    % make bounded
    X0 = Shapes.Sigma;
    X  = zclamp(X,0,Shapes.L0);
    
    k = 1;
    % construct shape-matrix for angular strains
    for ii = 1:Shapes.NModal(1)
        % not sure if interp1 is best/fastest option? maybe inverse lerp?
        Pc{k,1} = interp1(X0,Shapes.POD(:,ii),X);  
        k = k + 1;
    end
    
    % construct shape-matrix for linear strains
    for ii = 1:Shapes.NModal(2)
        Pc{k,1} = interp1(X0,Shapes.PODQ(:,ii),X);
        k = k + 1;
    end

    P = blkdiag(Pc{:});

end
%-------------------------------------------------- compute Cosserat string
function [Node, Rotation, Kappa, Gamma] = ReconstructField(Shapes,N,R,S)
    
P   = Shapes.Filter;  
lst = 1:Shapes.Fem.NNode;

Gamma = zeros(Shapes.NNode,1);
Kappa = zeros(Shapes.NNode,1);
Rotation = {};
Node     = P*N;

% loop over each node in curve to find geometric strains
for ii = 1:Shapes.NNode
    
    RRe = 0;
    UUe = 0;
    
    W   = P(ii,:);
    for jj = lst(abs(P(ii,:))>0)
        UUe = UUe + W(jj)*S{jj};
        RRe = RRe + W(jj)*R{jj};
    end
    
    [Ur,~,Vr] = svd(RRe);
    Re = (Ur*Vr.');
% 
%     Re = R{ii};
%     UUe = R{ii};
    
    Rotation{ii,1} = Re;
    Tangent(ii,:)  = Re(:,1).';
    Gamma(ii,1)    = sqrt(UUe(1,1));
end

%Kappa = [];

for ii = 2:Shapes.NNode-1
    Ni  = Node(ii,:);
    Nii = Node(ii+1,:);
    Nip = Node(ii-1,:);
    
    g  = [Rotation{ii,1},    [Ni(1);0;Ni(2)];   zeros(1,3), 1];
    dg = ([Rotation{ii+1,1}, [Nii(1);0;Nii(2)]; zeros(1,3), 1] ...
        - [Rotation{ii-1,1}, [Nip(1);0;Nip(2)]; zeros(1,3), 1])/(2*Shapes.ds);
    
    XI = inv(g)*dg;
    
    Kappa(ii,1) = -XI(1,2);
    %Gamma(ii,1) = XI(1,4);
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
%         Kappa(ii,1) = angle/(norm(Node(ii+1,:) - Node(ii,:)) + ...
%             norm(Node(ii,:) - Node(ii-1,:)));
%     end   
% 
% end


Kappa(1,1) = Kappa(2,1);
Gamma(1,1) = Gamma(2,1);
Kappa(end,1) = Kappa(end-1,1);
Gamma(end,1) = Gamma(end-1,1);
Kappa    = GaussianFilter(Kappa,round(Shapes.NNode/200));
Gamma    = GaussianFilter(Gamma,round(Shapes.NNode/200));

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
for ii = 1:Shapes.NNode
%     if ii < Shapes.NNode
%         dgam(ii,:) = Node(ii + 1,:) -  Node(ii,:);
%     else
%         dgam(ii,:) = Node(ii,:) - Node(ii - 1,:);
%     end
%     
    T(ii,:) = dgam(ii,:)/norm(dgam(ii,:));
    
end

I = null(round(T(1,:)));
Normal = I(:,1);
    
% compute curvature
for ii = 2:Shapes.NNode-1
    
    t1 = T(ii - 1,:);
    t2 = T(ii,:);
%     
    %t1 = (T(ii,:) + T(ii-1,:)); t1 = t1/norm(t1);
    %t2 = (T(ii,:) + T(ii+1,:)); t2 = t2/norm(t2);
    
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

R = Quat2Rot(Q);
A = StrainMap(R*Kap(:));

dg = zeros(7,1);

dg(1:4) = ((2*norm(Q))^(-1))*A*Q;
dg(5:7) = R*Gam(:);

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
rho = single(D(:,end)<1e-5);

I0 = reshape(rho,[N,N]);

% https://ocw.mit.edu/courses/aeronautics-and-astronautics/
% 16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec26.pdf
x0 = x0 - Shapes.Center(1);
y0 = y0 - Shapes.Center(2);
X0 = X0 - Shapes.Center(1);
Y0 = Y0 - Shapes.Center(2);

% evaluate slice volume
Shapes.Att   = sum(sum(I0*dv));

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
%----------------------------------------------------------- strain mapping
function Ktt = LinearStiffnessTensor(Shapes)
    E0 = Shapes.E;
    G0 = (E0)/(2*(1+Shapes.Nu));
  
    QE  = diag([E0,G0,G0]);
    QR  = diag([G0,E0,E0]);    
    Ktt = blkdiag(QR*Shapes.Jtt,QE*Shapes.Att);    
end

end
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