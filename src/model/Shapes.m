classdef Shapes

    properties (Access = public)
        Fem;
        BdBox;
        Table;
        Node;
        NNode;
        %NModal;
        %NDof;
        %NDim;
        %Phi;
        POD;
        PODEnergy;
        %Tracker;
        %posData;
        %xiData;
        %Sdomain;
    end
    
    properties (Access = private)
        Ba;
        xia0;
        Phi0;
        DiffStepSize;
        Filter;
        FilterRadius;
    end
   
%--------------------------------------------------------------------------
methods  
%----------------------------------------------- MODAL SHAPE RECONSTRUCTION
function obj = Shapes(Fem,Table,varargin) 
    
    obj.Table = Table; 
    obj.Fem   = Fem;  
    obj.xia0  = [0,0,0,1,0,0].';
    obj.BdBox = Fem.BdBox;
    obj.NNode = 300;
    obj.DiffStepSize = 1e-4;
    
    obj.FilterRadius = 0.02*mean([obj.BdBox(2)-obj.BdBox(1);...
                                  obj.BdBox(4)-obj.BdBox(3)]);
     
%     if ~isempty(Fem.get('Output'))
%         list = Fem.get('Output');
%         obj.Tracker = list(:,1);
%     else
%         error('Fem.Output is empty -- no tracker nodes selected in Fem');
%     end
%     
%     set = 1:6;
%     I6  = eye(6);
%     xa  = set(logical(Table));
   
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
%     
%     obj.Ba   = I6(:,xa);
%     obj.NDof = size(obj.Ba,2);
%     obj.Phi  = @(x) ShapeFunction(obj,x);
%     obj.NDim = obj.NDof*obj.NModal;
   
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
else, Request = varargin{1}; end

figure(101);

switch(Request)
    case('Base'),      Z = Shapes.POD;
    case('POD'),       Z = Shapes.POD;
    otherwise,         Z = Shapes.POD; Request = 'Base';
end
    
clf; axis equal; 
    
if strcmp(Request,'Base') || strcmp(Request,'POD') 
    subplot(2,1,1);
    if nargin>2,
        N = varargin{2};
    else
        N = 5;
    end
    
    X = linspace(0,1,size(Shapes.POD,1));
    
    for ii = 1:N
       plot(X,Shapes.POD(:,ii),'Color',col(ii),...
           'linewidth',1.5); hold on; grid on;
    end
    
    subplot(2,1,2);
    plot(X,(Shapes.PODEnergy).^(0.25),'-o','linewidth',1.5);
    grid on;
end

end
%---------------------------------------------------------------- show mesh
function Shapes = reference(Shapes,varargin)
X0 = varargin{1};
XL = varargin{2};

dX = XL(:) - X0(:);
if dX(2) < dX(1)
    Shapes.Node = [linspace(X0(1),XL(1),Shapes.NNode).',...
                   zeros(Shapes.NNode,1)];
else
    Shapes.Node = [zeros(Shapes.NNode,1),...
                   linspace(X0(2),XL(2),Shapes.NNode).'];
end

Shapes = GenerateRadialFilter(Shapes);

end
%---------------------------------------------------------------------- set
function Shapes = rebuild(Shapes,varargin)
    
    set = 1:6;
    I6 = eye(6);
    xa = set(logical(Shapes.Table));
   
    for ii = 1:2:length(varargin)
        Shapes.(varargin{ii}) = varargin{ii+1};
    end
    
    Shapes.Ba   = I6(:,xa);
    Shapes.NDof = size(Shapes.Ba,2);
    Shapes.Phi  = @(x) PODShapeFunction(Shapes,x);
    Shapes.NDim = Shapes.NDof*Shapes.NModal;
    
end     
%--------------------------------------------------------------------------
function Shapes = fit(Shapes)
fem = Shapes.Fem;
P = Shapes.Filter;
t = fem.Log.t;

for ii = 50:numel(t)
    N = fem.Log.Node{ii};
    
    Nfit = P*N;
    plot(Nfit(:,1),Nfit(:,2),'-x'); hold on;
    plot(fem.Log.Nx(:,ii),fem.Log.Ny(:,ii),'-o'); 
    axis equal;
 1   
end

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
function [p, g, X] = string(Shapes,q)
    
    g0 = [1,0,0,0,0,0,0];
    X = linspace(0,Shapes.Sdomain,Shapes.Quality);
    
    [~,y] = ode23(@(t,x) ForwardODE(Shapes,t,x,q),X,g0);

    g = y(:,1:7);
    p = [y(:,7),y(:,6),y(:,5)];
    
end
end

methods (Access = private)
%-------------------------------------------------- compute Cosserat string
function Shapes = GenerateRadialFilter(Shapes)

PS2 = Shapes.Node;    
PS1 = Shapes.Fem.get('Node0');

d = DistancePointSet(PS1,PS2,Shapes.FilterRadius);

if ~isempty(d)
    P = sparse(d(:,1),d(:,2),1-d(:,3)/Shapes.FilterRadius,...
        Shapes.NNode,Shapes.Fem.NNode);
    P = spdiags(1./sum(P,2),0,size(P,1),size(P,1))*P;
else
    P = 1;
end

Shapes.Filter = P;

end
%---------------------------------------------------------------------- set
function P = ShapeFunction(Shapes,X)

    P11 = zeros(1,max(Shapes.NModal));
    Pc  = cell(Shapes.NDof,1);
    X   = single(X)/Shapes.Sdomain;

    for ii = 1:max(Shapes.NModal)
        P11(1,ii) = legendre(X,ii-1);
    end

    P = P11;

    for ii = 1:Shapes.NDof
        Pc{ii,1} = P;
    end

    P = blkdiag(Pc{:})+ 1e-26*X;

end
%---------------------------------------------------------------------- set
function P = PODShapeFunction(Shapes,X)
    
    P11 = zeros(1,max(Shapes.NModal));
    Pc  = cell(Shapes.NDof,1);
    X   = single(X)/Shapes.Sdomain;
    N   = length(Shapes.POD(:,1));

    for ii = 1:max(Shapes.NModal)
        id0 = clamp(floor(X*N),1,80);
        id1 = clamp(ceil(X*N),1,80);
        p = Shapes.POD(:,ii);
        
        P11(1,ii) = lerp(p(id0),p(id1),((X*N)-...
            floor(X*N)));
    end

    P = P11;

    for ii = 1:Shapes.NDof
        Pc{ii,1} = P;
    end

    P = blkdiag(Pc{:})+ 1e-26*X;
end
%--------------------------------------- forwards integration of kinematics
function dg = ForwardODE(Shapes,t,g,q)
ee = Shapes.Ba*Shapes.Phi(t)*q + Shapes.xia0(:);

Kappa = ee(1:3);
Gamma = ee(4:6);
Q     = g(1:4);

R = Quat2Rot(Q);
A = StrainMap(R*Kappa(:));

dg = zeros(7,1);

dg(1:4)   = ((2*norm(Q))^(-1))*A*Q;
dg(5:7)   = R*Gamma(:);

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
%---------------------------------------------------------- material filter 
function d = DistancePointSet(PS1,PS2,R)

d = cell(size(PS1,1),1);
    
for el = 1:size(PS1,1)   
    
    if size(PS1,2) == 3
        dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2 + ...
            (PS1(el,3)-PS2(:,3)).^2 );    
    else    
        dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^6);
    end

    [I,J] = find(dist<=R);   
    d{el} = [I,J+(el-1),dist(I)];
end

% matrix of indices and distance value
d = cell2mat(d); 
end



