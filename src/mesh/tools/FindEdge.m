function x = FindEdge(Mesh,varargin)

tol   = BuildTolerance(Mesh.Node); 
BdBox = BoundingBox(Mesh.Node);

Request = varargin{1}; 

switch(Request)
    case('Hole');       x = FindHoles(Mesh,varargin{end});
    case('TopHole');    x = FindTopHoles(Mesh,tol);
    case('EdgeSelect');   x = FindEdgeSelect(Mesh,BdBox,varargin{2:end});
end

end

function id = FindHoles(Mesh,List)
Bnd = Mesh.get('BndMat');
for ii = 1:length(List)
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

function id = FindEdgeSelect(Mesh,BdBox,Point,AngleMax)

if nargin < 4, AngleMax = 90 - 1e-6; end
Bnd = Mesh.get('BndMat');

for ii = 1:length(Bnd)
    E = unique(Bnd{ii}(:),'stable');
    E = [E;E(1)];
    S{ii,1} = E;
end   

id = Location(Mesh.Node,BdBox,Point);

if length(S) == 1, 
   Loop = S{1};
else
   for ii = 1:length(S);
        if ~ismember(S{ii},id)
        else, break;
        end
   end
   Loop = S{ii};
end

V = Mesh.Node(Loop,:);
dsx = diff(V(:,1));
dsy = diff(V(:,2));
ds = sqrt(dsx.^2+dsy.^2);
Tx = dsx./ds;
Ty = dsy./ds;

% Second derivative & curvature
ds2 = 0.5*(ds([end,1:end-1])+ds);
Hx = diff(Tx([end,1:end]))./ds2;
Hy = diff(Ty([end,1:end]))./ds2;
x = V(1:end-1,1);
y = V(1:end-1,2); % remove repeated point

Angles = sqrt(Hx.^2+Hy.^2).*ds*(180/pi);

index0 = find(Loop == id);
n = 1; m = 1;

%forward loop
for ii = (index0+1):1:length(Loop)
    if Angles(ii) <= AngleMax
        n = n+1;
    else, 
        break 
    end
end
%backwards loop
for ii = (index0-1):-1:1
    if Angles(ii) <= AngleMax
        m = m+1;
    else, 
        break; 
    end
end

id = {flipud(Loop((index0-m):1:(index0+n)))};

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

