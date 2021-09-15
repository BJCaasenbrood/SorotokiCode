% adaption of the code: https://stackoverflow.com/questions/6485908/basic-dual-contouring-theory
function [F,V] = DualContour(Sdf,N)
if nargin < 2
   N = 16; 
end
tic;
F = []; 
V = [];

deps = 1e-3;

% set Sdfs
f  = @(x) Sdf.eval(x);
df = @(x) gradientEstimateMD(f,x,deps); % aka normal estimation

% generate unit cell
dx = (Sdf.BdBox(2) - Sdf.BdBox(1))/(N-1);
dy = (Sdf.BdBox(4) - Sdf.BdBox(3))/(N-1);
dz = (Sdf.BdBox(6) - Sdf.BdBox(5))/(N-1);

[Uv, Ue] = unitCube(dx,dy,dz);

% generate meshgrid of cube centers
[X,Y,Z] = meshgrid(linspace(Sdf.BdBox(1),Sdf.BdBox(2),N),...
                   linspace(Sdf.BdBox(3),Sdf.BdBox(4),N),...
                   linspace(Sdf.BdBox(5),Sdf.BdBox(6),N));

% generate [8x3xN^3] matrix of cube vertices
Vmat = repmat(Uv,1,1,N^3) + ...
    repmat(reshape([X(:),Y(:),Z(:)].',1,3,N^3),8,1,1);

%% check for sign changes
Vopt    = zeros(N^3,3);
Id      = zeros(N^3,1);
Sdetect = false(N^3,1);
Edetect = false(12,N^3);

% opt     = optimoptions('lsqlin','Display','none',...
%     'MaxIterations',1e3,'Algorithm','trust-region-reflective',...
%     'ConstraintTolerance',1e-8,...
%     'OptimalityTolerance',1e-8); %'');

opt = optimoptions(@lsqlin,'Algorithm','active-set',...
    'MaxIterations',1,'Display','off');

for ii = 1:size(Vmat,3)
    D = f(Vmat(:,:,ii)); 
    S = sign(D(:,end));
    
    if abs((abs(S)).'*S) < 8
        Sdetect(ii) = true;
        De  = D(Ue(:),end); 
        Ve  = Vmat(Ue(:),:,ii);
        V1  = Ve(1:2:24,:); V2 = Ve(2:2:24,:);  % location edge vertices
        dV1 = De(1:2:24); dV2  = De(2:2:24);    % sdf at edge vertices
        Se = logical(sign(dV1) - sign(dV2));    % find edges with sign change
        
        %backwards and forward lerp
        t = -dV1(Se)./(dV2(Se)   - dV1(Se));    % backwards
        P = V1(Se,:) + (V2(Se,:) -V1(Se,:)).*t; % forwards: intersec. points
        
        Gi = df(P);                % find gradients
        Ni = Gi./vecnorm(Gi.').';  % normalize vectors
        
        % solve Quadradic Error Function (QEF)
        A = Ni;
        b = dot(Ni,P,2); 
        lb = [min(Ve(:,1)); min(Ve(:,2)); min(Ve(:,3))];
        ub = [max(Ve(:,1)); max(Ve(:,2)); max(Ve(:,3))];
        
        %x0 = ;
        y = lsqlin(A,b,[],[],[],[],lb,ub,mean(P,1),opt);
        %y = A\b;
        Vopt(ii,:) = (y).';   % optimal vertex placement in unit-cube
        Id(ii,1)   = ii;
        Edetect(:,ii) = Se;
    end
end

Clist = 1:(N^3); 
Clist = Clist(Sdetect);
Elist = 1:12;

% plot3(Vopt(:,1),Vopt(:,2),Vopt(:,3),'b.');
% axis equal
% grid on;
% hold all;
% view(30,30);

% SEMI-IDEA how ;)
for ii = Clist
    edge = Elist(Edetect(:,ii));
    
    for jj = edge

        [yi,xi,zi] = ind2sub([N,N,N],ii);
        Cids = edgeToCubeLoop(jj) + [xi,yi,zi];
        
        Ftmp = zeros(1,4);
        
        for kk = 1:4
             if (Cids(kk,1) ~= 0) && (Cids(kk,2) ~= 0) && (Cids(kk,3) ~= 0)
                Ftmp(kk) = sub2ind([N,N,N],Cids(kk,2),...
                                           Cids(kk,1),...
                                           Cids(kk,3));
             end
        end
        
        
        % construct quad;
        if length(unique(Ftmp)) == 4        
            %F = [F;Ftmp(1),Ftmp(2),Ftmp(3);
            %       Ftmp(1),Ftmp(3),Ftmp(4)];
            F = [F;Ftmp];
        end
      
    end
    
end

%Vopt(Id == 0,:) = NaN;
[V,~,Ic] = unique(Vopt,'rows');
F = Ic(F);
%F = unique(sort(F,2),'rows');
%V = Vopt;

end
%--------------------------------------------- estimate the gradient of SDF
function df = gradientEstimateMD(f,x,delta)
df = zeros(length(x),3);

for ii = 1:3
    xmin = x; xmax = x;
    xmin(:,ii) = x(:,ii) - delta;
    xmax(:,ii) = x(:,ii) + delta;
    Dmin = f(xmin); Dmin = Dmin(:,end);
    Dmax = f(xmax); Dmax = Dmax(:,end);
    df(:,ii) = (Dmax - Dmin)/2/delta;
end

end
%--------------------------------- generate vertices and edges of unit-cube
function [V,E] = unitCube(dx,dy,dz)
V = unique(nchoosek([0,1,0,1,0,1],3),'rows');
E = unique(nchoosek(1:8,2),'rows');

% elimate all diagonal edges
X     = V(E(:,1),:) - V(E(:,2),:);
edgeL = diag(X*X.');
E = E(edgeL == 1,:);

V = V.*[dx,dy,dz];
E = E';
end

%----------------------------------------------------------------------
function lp = edgeToCubeLoop(id)

X = [1,0,0];
Y = [0,1,0];
Z = [0,0,1];
O = [0,0,0];

switch id
    case 1,  lp = [O;-X;-X-Y;-Y];
    case 9,  lp = [O;-Y;+X-Y;+X];
    case 12, lp = [O;+X;+X+Y;+Y];
    case 6,  lp = [O;+Y;-X+Y;-X]; %_
    case 2,  lp = [O;-Z;-X-Z;-X];
    case 4,  lp = [O;-X;-X+Z;+Z];
    case 11, lp = [O;+Z;+X+Z;+X];
    case 10, lp = [O;+X;+X-Z;-Z];  
    case 3,  lp = [O;-Y;-Y-Z;-Z];
    case 7,  lp = [O;-Z;+Y-Z;+Y];
    case 8,  lp = [O;+Y;+Y+Z;+Z];
    case 5,  lp = [O;+Z;-Y+Z;-Y];
    otherwise, lp = [O;O;O;O];
end

end

%----------------------------------------------------------------------
function F = quad2tri(F0)

nel1 = size(F0,1) ;        
nel2 = 2*nel1 ;
F0 = zeros(nel2,3) ;
%
F0(1:2:end,:) = F0(:,1:3);
F0(2:2:end,:) = F0(:,[1 3 4]);

end