classdef Mesh

    properties (Access = public)
        SDF;
        BdBox;
        Node;
        Element;
        NNode;
        NElem;
        Center;
        Dimension = 2;
    end
    
    properties (Access = private)
        Area;
        Support;
        Load;
        Spring;
        Output;
        OutputDir;
        PressureCell;
        FixedDensity;
        Convergence;
        Velocity;
        Adjecency;
        FaceToNode;
        NodeToFace;
        Boundary;
        Normal;
        Iteration = 0;
        ElemMat = -1;
        eps = 1e-5; 
        eta = 0.9; 
        Triangulate = false;
        MaxIteration = 100;
        ShowMeshing = false;
        CollapseTol = 0.1;
    end
    
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = Mesh(SDF,varargin) 
    obj.SDF = SDF;
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- get     
function varargout = get(Mesh,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Mesh.(varargin{ii});
        end
    else
        varargout = Mesh.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function Mesh = set(Mesh,varargin)
    for ii = 1:2:length(varargin)
        Mesh.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- set
function Mesh = generateMesh(Mesh)
Mesh.Iteration = 1; 
Mesh.Convergence = zeros(Mesh.MaxIteration,1);

if isempty(Mesh.Center)
    Pc = randomPointSet(Mesh); 
else
    Mesh.MaxIteration = 2;
    Pc = Mesh.Center; 
    d = Mesh.SDF(Pc);
    Pc = Pc(d(:,end)<0,:);  
    Mesh.NElem = length(Pc);
end

A_ = (Mesh.BdBox(2)-Mesh.BdBox(1))*...
     (Mesh.BdBox(4)-Mesh.BdBox(3));

while Mesh.Iteration < Mesh.MaxIteration
  
  P = Pc;
  
  Rp = pointSetReflect(Mesh,P,A_); 
  
  Rb = boundingReflect(Mesh);

  [v,f] = voronoin([P;Rp;Rb],{'Qt','Pp'});   
  
  [Pc,A] = computeCentroid(Mesh,f,v); 
  
  Mesh.Velocity = vecnorm((P - Pc)')';
  
  A_ = sum(abs(A));
  
  Mesh.Convergence(Mesh.Iteration) = sqrt(sum((A.^2).*sum((Pc-P)...
      .*(Pc-P),2)))*Mesh.NElem/(A_^1.5);
  
  Mesh.Iteration = Mesh.Iteration + 1;
  
  if Mesh.ShowMeshing
     Mesh.Center = Pc;
     Mesh.Node = v;
     Mesh.Element = f(1:Mesh.NElem);
     show(Mesh,'Velocity');
  end
end

f = f(1:Mesh.NElem);

[v,f] = ExtractNode(Mesh,v,f);
[v,f] = CollapseEdges(Mesh,v,f);
[v,f] = ResequenceNodes(Mesh,v,f);

if Mesh.Triangulate
    Mesh.Node = v;
    Mesh.Element = f;
    Mesh.NNode = length(v);
    Mesh.NElem = length(f);
    [v,f] = MeshTriangulation(Mesh,P,v,f);
    [v,f] = ResequenceNodes(Mesh,v,f);
    Mesh.NElem = length(f);
    Mesh.NNode = length(v);
end

[Pc,A] = computeCentroid(Mesh,f,v); 
normal = computeNormalVector(Mesh,f,v,Pc);

Mesh.Center = Pc;
Mesh.Node = v;
Mesh.Element = f;
Mesh.NNode = length(v);
Mesh.NElem = length(f);
Mesh.Area = A;
Mesh.Normal = normal;

Mesh = ElementAdjecency(Mesh);
end

%---------------------------------------------------------- add constraints
function Mesh = addConstraint(Mesh,varargin)
    
for ii = 1:2:length(varargin)
  if size(varargin{ii+1},2) == 3
      Mesh.(varargin{ii}) = varargin{ii+1};
  else
      warning([varargin{ii}, ' has incorrect input'] );
  end
end
    
end

%---------------------------------------------------------------- show mesh
function show(Mesh,varargin)
    
if nargin<2, Request = -1; 
else, Request = varargin{1}; end

Mesh = ElementAdjecency(Mesh);

switch(Request)
    case('Velocity'), Z = (Mesh.Velocity')';  caxis([0 0.01]);
    case('Gradient'), Z = GradientField(Mesh,Mesh.Node);
    case('Node'), Z = varargin{2};
    case('Element'), Z = varargin{2};
    otherwise; Dist = Mesh.SDF(Mesh.Node); Z = Dist(:,end);
end

if ~strcmp(Request,'Node') && ~strcmp(Request,'Element')
clf; axis equal; axis off; hold on;
    
patch('Faces',Mesh.ElemMat,'Vertices',Mesh.Node,...
    'FaceVertexCData',Z,'Facecolor','flat','LineStyle','-',...
    'Linewidth',1.0,'FaceAlpha',1.0,'EdgeColor','k');

patch('Faces',Mesh.Boundary,'Vertices',Mesh.Node,...
    'LineStyle','-','Linewidth',2,'EdgeColor','k');
elseif strcmp(Request,'Node')
hold on;
plot(Mesh.Node(Z,1),Mesh.Node(Z,2),'.','Color','r');
elseif strcmp(Request,'Element')
hold on;
plot(Mesh.Center(Z,1),Mesh.Center(Z,2),'.','Color','r');    
end

axis(Mesh.BdBox);
axis tight;
drawnow;
end

%--------------------------------------------------------------- find nodes
function NodeList = FindNodes(Mesh,Request)
    n = Mesh.Node;
    NodeList = FindNodes(n,Request);
end

%--------------------------------------------------------------- find nodes
function ElementList = FindElements(Mesh,Request)
    n = Mesh.Center;
    ElementList = FindNodes(n,Request);
end
%-------------------------------------------------------------- END METHODS
end

methods (Access = private)
%----------------------------------------------- generate a random pointset
function P = randomPointSet(Mesh)
P = zeros(Mesh.NElem,Mesh.Dimension);
B = Mesh.BdBox;
Ctr = 0;
while(Ctr < Mesh.NElem)
    Y = zeros(Mesh.NElem,Mesh.Dimension);
    Y(:,1) = (B(2)-B(1))*rand(Mesh.NElem,1)+B(1);
    Y(:,2) = (B(4)-B(3))*rand(Mesh.NElem,1)+B(3);
    d = Mesh.SDF(Y);
    I = find(d(:,end)<0);               
    NumAdded = min(Mesh.NElem-Ctr,length(I));
    P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
    Ctr = Ctr+NumAdded;
end
end 

%--------------------------------------------------------- reflect pointset
function Rp = pointSetReflect(Mesh,P,A)
Alpha = 1.5*sqrt(A/Mesh.NElem);
d = Mesh.SDF(P);  

% number of constituent bdry segments
NBdrySegs = size(d,2)-1;       
n1 = (Mesh.SDF(P+repmat([Mesh.eps,0],Mesh.NElem,1))-d)/Mesh.eps;
n2 = (Mesh.SDF(P+repmat([0,Mesh.eps],Mesh.NElem,1))-d)/Mesh.eps;

%Logical index of seeds near the boundary
I  = abs(d(:,1:NBdrySegs))<Alpha; 

P1 = repmat(P(:,1),1,NBdrySegs); 
P2 = repmat(P(:,2),1,NBdrySegs); 
Rp(:,1) = P1(I)-2*n1(I).*d(I);  
Rp(:,2) = P2(I)-2*n2(I).*d(I);

d_R_P = Mesh.SDF(Rp);
J    = abs(d_R_P(:,end))>=Mesh.eta*abs(d(I)) & d_R_P(:,end)>0;
Rp   = Rp(J,:); 
Rp   = unique(Rp,'rows');
end

%------------------------------------------------ compute centroid polygons
function Rb = boundingReflect(Mesh)
    
tmp = Mesh.BdBox; 

a = 0.5;
tmp(1) = Mesh.BdBox(1)-a*( Mesh.BdBox(2) - Mesh.BdBox(1)); 
tmp(2) = Mesh.BdBox(2)+a*( Mesh.BdBox(2) - Mesh.BdBox(1)); 
tmp(3) = Mesh.BdBox(3)-a*( Mesh.BdBox(4) - Mesh.BdBox(3)); 
tmp(4) = Mesh.BdBox(4)+a*( Mesh.BdBox(4) - Mesh.BdBox(3)); 

x = linspace(tmp(1),tmp(2),2);
y = linspace(tmp(3),tmp(4),2);

[X,Y] = meshgrid(x,y);
Rb = [X(:),Y(:)];  

end

%------------------------------------------------ compute centroid polygons
function [Pc,A] = computeCentroid(Mesh,f,v)

Pc = zeros(Mesh.NElem,2); 
A  = zeros(Mesh.NElem,1);

for el = 1:Mesh.NElem

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

%--------------------------------------------- compute normal from centroid
function n = computeNormalVector(Mesh,f,v,CNode)

N = Mesh.NElem;
n{N} = [];

for i = 1:length(f)
    id = f{i};
    Nd = v(id,:) - repmat(CNode(i,:),length(id),1);
    n{i} = -Nd;
end

end

%----------------------------------------------- triangulate polygonal mesh
function [v, f] = MeshTriangulation(Mesh,Center,v0,f0)

f = [];

for ii = 1:Mesh.NElem
    el = f0{ii};
    n = numel(el);
    elem = [el(1:n)', [el(2:n)'; el(1)], ...
            repmat(Mesh.NNode+ii,n,1)];
    f = [f; elem];
end

f = num2cell(f,2);
v = [v0;Center];
end

%--------------------------------------------------- collapse smaller edges
function [Node0,Element0] = CollapseEdges(Mesh,Node0,Element0)

while(true)
  cEdge = [];
  for el=1:size(Element0,1)
    %Cannot collapse triangles
    if size(Element0{el},2)<4, continue; end

    vx = Node0(Element0{el},1); 
    vy = Node0(Element0{el},2); 
    nv=length(vx);
    beta = atan2(vy-sum(vy)/nv, vx-sum(vx)/nv);
    beta = mod(beta([2:end 1]) -beta,2*pi);
    betaIdeal = 2*pi/size(Element0{el},2);
    Edge = [Element0{el}',Element0{el}([2:end 1])'];
    cEdge = [cEdge; Edge(beta<Mesh.CollapseTol*betaIdeal,:)];
  end

  if (size(cEdge,1)==0)
    break; 
  end

  cEdge = unique(sort(cEdge,2),'rows');
  cNode = 1:size(Node0,1);

  for i=1:size(cEdge,1)
    cNode(cEdge(i,2)) = cNode(cEdge(i,1));
  end

  [Node0,Element0] = Rebuild(Mesh,Node0,Element0,cNode);
end

end

%------------------------------------------------------------------ rebuild
function [Node,Element] = Rebuild(Mesh,Node0,Element0,cNode)

Element   = cell(Mesh.NElem,1);
[~,ix,jx] = unique(cNode);

if ~isequal(size(jx),size(cNode)), jx=jx'; end 
if size(Node0,1)>length(ix), ix(end)=max(cNode); end

Node = Node0(ix,:); 

for el=1:size(Element0,1)
  Element{el} = unique(jx(Element0{el}));
  vx = Node(Element{el},1); 
  vy = Node(Element{el},2); 
  nv = length(vx);

  [~,iix] = sort(atan2(vy-sum(vy)/nv,vx-sum(vx)/nv));
  Element{el} = Element{el}(iix);
end

end

%------------------------------------------------------- extract nodal data
function [Node,Element] = ExtractNode(Mesh,Node0,Element0)
map = unique([Element0{1:Mesh.NElem}]);
cNode = 1:size(Node0,1);
tmp = setdiff(cNode,map);
cNode(tmp) = max(map);

[Node,Element] = Rebuild(Mesh,Node0,Element0(1:Mesh.NElem),cNode);
end

%---------------------------------------------------------- resequence mesh
function [Node,Element] = ResequenceNodes(Mesh,Node0,Element0)
NNode0=size(Node0,1); NElem0=size(Element0,1);
ElemLnght = cellfun(@length,Element0); 

nn = sum(ElemLnght.^2); 
i = zeros(nn,1); 
j = zeros(nn,1); 
s = zeros(nn,1); 

index=0;

for el = 1:NElem0
  eNode=Element0{el}; ElemSet=index+1:index+ElemLnght(el)^2;
  i(ElemSet) = kron(eNode,ones(ElemLnght(el),1))';
  j(ElemSet) = kron(eNode,ones(1,ElemLnght(el)))';
  s(ElemSet) = 1;
  index = index+ElemLnght(el)^2;
end

K = sparse(i,j,s,NNode0, NNode0);
p = symrcm(K);
cNode(p(1:NNode0))=1:NNode0;
[Node,Element] = Rebuild(Mesh,Node0,Element0,cNode);

end

%------------------------------------------------ generate elemental matrix
function ElemMat = GenerateElementalMatrix(Mesh)
El = Mesh.Element(1:Mesh.NElem)';                 
MaxNVer = max(cellfun(@numel,El));      
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))]; 
ElemMat = cellfun(PadWNaN,El,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});       
end

%----------------------------------------------- generate elemental adjency
function Mesh = ElementAdjecency(Mesh)
face = Mesh.Element;
edges = cellfun(@numel,face);
n = Mesh.NElem;    
p = max(cellfun(@max,face));    
MaxNVer = max(cellfun(@numel,face));     
tri = GenerateElementalMatrix(Mesh);
set = 1:n;

maxver = num2cell(1:size(tri,2));
fnc = @(x) sum(tri(:,x)>=1);
countsperrow = cellfun(@(x) fnc(x), maxver);

idx = 0;
I = zeros(sum(countsperrow),1);
J = zeros(sum(countsperrow),1);

for ii = 1:MaxNVer
   id = ~isnan(tri(:,ii));
   i = set(id)';
   j = tri(id,ii);
   I(idx+1:idx+countsperrow(ii)) = i;
   J(idx+1:idx+countsperrow(ii)) = j;
   idx = idx + countsperrow(ii);
end

M = sparse(I,J,1,n,p);

C = M*M'-diag(diag(M*M'));
EdgeMat = zeros(sum(edges),2);
idx = 0;

for ii = 1:Mesh.NElem
    [tmp,N]= EdgeMatrix(Mesh,face{ii});
    EdgeMat(idx+1:N+idx,:) = [face{ii}(tmp(:,1))',...
                              face{ii}(tmp(:,2))'];
    idx = idx + N;
end

E = sort(EdgeMat')';
[u,~,n] = unique(E,'rows');
counts = accumarray(n(:), 1);
B = u(counts==1,:);

Mesh.Adjecency = sparse(double(C>1));
Mesh.NodeToFace = (M./sum(M,1))';
Mesh.FaceToNode  = ((M')./sum(M,2)')';
Mesh.Boundary = B;
Mesh.ElemMat = tri;
end

%------------------------------------------------------ compute edge matrix
function [e,n] = EdgeMatrix(Mesh,face)
q = Mesh.NElem;
n = numel(face);
v = 1:n;
e = [v(:),[v(2:end)';v(1)]];
end

%------------------------------------------------------- get gradient field
function RGB = GradientField(Mesh,P)
    
d = Mesh.SDF(P);  
N = length(P);
n1 = (Mesh.SDF(P+repmat([Mesh.eps,0],N,1))-d)/Mesh.eps;
n2 = (Mesh.SDF(P+repmat([0,Mesh.eps],N,1))-d)/Mesh.eps; 

RGB = atan2(n2(:,end),n1(:,end));
end

end
end

