function x = FindEdge(Mesh,varargin)
%FINDEDGE Returns a cell-list of nodes that form an edge.
%
%   x = FINDEDGE(Mesh, 'Hole', A)   -  Returns holes of index A
%   x = FINDEDGE(Mesh, 'TopHole')   -  Returns holes at the top
%   x = FINDEDGE(Mesh, 'AllHole')   -  Returns all holes
%
%   x = FINDEDGE(Mesh, 'EdgeSelect', P0, alpha) -  Returns an edge closest
%       to point P0, and neighbouring edges within an angle limit <= alpha.
%       alpha by default 90 degrees.
%
%   See also FINDNODE
%   Brandon Caasenbrood
%   2020, MIT LICENSE.

tol   = BuildTolerance(Mesh.Node); 
BdBox = BoundingBox(Mesh.Node);

Request = varargin{1}; 

switch(Request)
    case('Hole');       x = FindHoles(Mesh,varargin{end});
    case('TopHole');    x = FindTopHoles(Mesh,tol);
    case('AllHole');    x = FindAllHoles(Mesh);
    case('BoxHole');    x = BoxHoles(Mesh,tol,varargin{2:end});
    case('EdgeSelect'); x = FindEdgeSelect(Mesh,BdBox,varargin{2:end});
end

end

function id = FindHoles(Mesh,List)
if ischar(List)
    id = FindAllHoles(Mesh);
    return
end

Bnd = Mesh.get('BndMat');
for ii = 1:numel(List)
    E = unique(Bnd{List(ii)}(:),'stable');
    E = [E;E(1)];
    id{ii,1} = E;
end   
end

function id = FindTopHoles(Mesh,tol)
Bnd = Mesh.get('BndMat');
Nds = Mesh.Node;

for ii = 1:length(Bnd)
    E = unique(Bnd{ii}(:),'stable');
    E = [E;E(1)];
    S{ii,1} = E;
end   

C = PolygonCenter(Nds,S);
id = S(abs(C(:,2)-max(C(:,2))) < tol);
end

function id = BoxHoles(Mesh,tol,Line)
Bnd = Mesh.get('BndMat');
Nds = Mesh.Node;

for ii = 1:length(Bnd)
    E = unique(Bnd{ii}(:),'stable');
    E = [E;E(1)];
    S{ii,1} = E;
end   

C = PolygonCenter(Nds,S);
d = dRectangle(C,Line(1)-eps,Line(2)+eps,Line(3)-eps,Line(4)+eps);
id = S(d(:,end)<tol);
end


function id = FindAllHoles(Mesh)
Bnd = Mesh.get('BndMat');
Nds = Mesh.Node;

for ii = 1:length(Bnd)
    E = unique(Bnd{ii}(:),'stable');
    E = [E;E(1)];
    S{ii,1} = E;
end   

id = S(2:end);
end

function id = FindEdgeSelect(Mesh,BdBox,Point,AngleMax)

if nargin < 4, AngleMax = 90 - 1e-6; end
Bnd = Mesh.get('BndMat');

for ii = 1:length(Bnd)
    E = unique(Bnd{ii}(:),'stable');
    E = [E;E(1)];
    S{ii,1} = E;
end   

id = Location(Mesh.Node,BdBox,Point);

if length(S) == 1
   Loop = S{1};
else
   for ii = 1:length(S)
        if ~ismember(S{ii},id)
        else, break;
        end
   end
   Loop = S{ii};
end

V = Mesh.Node(Loop,:);
% dsx = diff(V(:,1));
% dsy = diff(V(:,2));
% ds = sqrt(dsx.^2+dsy.^2);
% Tx = dsx./ds;
% Ty = dsy./ds;
% 
% % Second derivative & curvature
% ds2 = 0.5*(ds([end,1:end-1])+ds);
% Hx = diff(Tx([end,1:end]))./ds2;
% Hy = diff(Ty([end,1:end]))./ds2;
% x = V(1:end-1,1);
% y = V(1:end-1,2); % remove repeated point
% 
% Angles = sqrt(Hx.^2+Hy.^2).*ds*(180/pi);
% 
Angles = abs(DifferentialGeometry(V));

index0 = find(Loop == id);
n = 1; m = 1;

%forward loop
for ii = (index0+1):1:length(Loop)-1
    if Angles(ii) <= AngleMax
        n = n+1;
    else
        break 
    end
end

%backwards loop
for ii = (index0-1):-1:1
    if Angles(ii) <= AngleMax
        m = m+1;
    else
        break; 
    end
end

id = {flipud(Loop((index0-m+1):1:(index0+n)))};

end

function id = Location(Node,BdBox,Node0)
if length(BdBox) == 4
    Loc = sqrt((Node(:,1)-Node0(1)).^2+(Node(:,2)-Node0(2)).^2);
elseif length(BdBox) == 6
    Loc = sqrt((Node(:,1)-Node0(1)).^2+ ...
        (Node(:,2)-Node0(2)).^2 + (Node(:,3)-Node0(3)).^2);
end

[~,Loc] = sort(Loc); 
id = Loc(1);
end

function eps = BuildTolerance(Node)
BdBox = BoundingBox(Node);
eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
end

function BdBox = BoundingBox(Node)
if size(Node,2) == 2
BdBox = zeros(4,1);
BdBox(1) = min(Node(:,1));
BdBox(2) = max(Node(:,1));
BdBox(3) = min(Node(:,2));
BdBox(4) = max(Node(:,2));
else
BdBox = zeros(6,1);
BdBox(1) = min(Node(:,1));
BdBox(2) = max(Node(:,1));
BdBox(3) = min(Node(:,2));
BdBox(4) = max(Node(:,2));   
BdBox(5) = min(Node(:,3));
BdBox(6) = max(Node(:,3));  
end
end

function [Pc,A] = PolygonCenter(v,f)
N  = length(f);
Pc = zeros(N,2); 
A  = zeros(N,1);

for el = 1:length(f);

  vx = v(f{el},1); 
  vy = v(f{el},2); 
  nv = length(f{el}); 
  
  vxS = vx([2:nv 1]); 
  vyS = vy([2:nv 1]); 
  
  tmp = vx.*vyS - vy.*vxS;
  A(el) = 0.5*sum(tmp);
  Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*tmp),...
                            sum((vy+vyS).*tmp)];
end
end

function Theta = DifferentialGeometry(Node)
% http://page.math.tu-berlin.de/~bobenko/Lehre/Skripte/DDG_Lectures.pdf
N = length(Node);
T = zeros(N,3); 
Node0 = Node;
%dgam = zeros(N,3); 
%phi  = zeros(N,1);

[~, Fy] = gradient(Node);   
dgam = [Fy(:,1),Fy(:,2)]; 

[~, Fy] = gradient(Node0);   
dgam0 = [Fy(:,1),Fy(:,2)]; 

dl  = sqrt(sum(dgam.^2,2));
dl0 = sqrt(sum(dgam0.^2,2));
ds  = mean(dl0);

Gamma = dl./dl0;

% compute tangents
for ii = 1:N 
    T(ii,[1,3]) = dgam(ii,:)/norm(dgam(ii,:));
end

I = null(round(T(1,:)));
Normal = I(:,1);
    
% compute curvature
for ii = 2:N-1
    
    t1 = T(ii - 1,:);
    t2 = T(ii,:);
%     
    dir = sign(dot(so3(t2)*t1(:),Normal));
    angle = real(2*dir*acos(dot(t2,t1)));
    
    
    Kappa(ii,1) = angle/(norm(Node(ii+1,:) - Node(ii,:)) + ...
        norm(Node(ii,:) - Node(ii-1,:)));
end

Kappa(1,:) = Kappa(2,:) + ds*(Kappa(3,:) - Kappa(2,:));
Kappa(N,:) = Kappa(end-1,:) + ds*(Kappa(end,:) - Kappa(end-1,:));
Theta = (Kappa*ds)*(180/pi);

end
