classdef Mesh

    properties (Access = public)
        SDF;
        BdBox;
        Node;
        Element;
        NNode;
        NElem;
        Dim;
        Type;
        Center;
    end
    
    properties (Access = private)
        Area;
        Convergence;
        Velocity;
        Adjecency;
        FaceToNode;
        NodeToFace;
        ElementToFace;
        Boundary;
        Normal;
        Edge;
        Iteration;
        ElemMat;
        eps; 
        eta; 
        Triangulate;
        MaxIteration;
        Movie;
        MovieStart;
        CollapseTol;
        ConvNorm;
        ShowProcess;
        Colormap;
        LineStyle;
    end
    
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = Mesh(Input,varargin) 
    
    if isa(Input,'char')
       [~,~,ext] = fileparts(Input);
       if strcmp(ext,'.stl'), [f,v] = stlreader(Input);
       elseif strcmp(ext,'.obj'), [v,f] = objreader(Input);
       else, cout('err','* extension not recognized');
       end
       
       obj.Node    = v;
       obj.NNode   = length(v);
       obj.Element = num2cell(f,2);
       obj.NElem   = length(f);
       obj.BdBox   = boxhull(v);
       
    elseif isa(Input,'function_handle')
       obj.SDF = Input;
       obj.NElem = 200;
    end
    
    obj.ShowProcess  = false;
    obj.MaxIteration = 100;
    obj.Iteration    = 0;
    obj.ElemMat      = -1;
    obj.eps          = 1e-6; 
    obj.eta          = 0.9; %was 0.9
    obj.Triangulate  = false;
    obj.Movie        = false;
    obj.MovieStart   = false;
    obj.CollapseTol  = 0.2;
    obj.ConvNorm     = 1e-3;
    obj.ShowProcess  = false;
    obj.Colormap     = turbo;
    obj.LineStyle    = '-';
    obj.Type         = 'C2PX';
    
    for ii = 1:2:length(varargin)
        if strcmp(varargin{ii},'Quads')
            N = num2cell(varargin{ii+1});
            obj.Center = Quads(obj.BdBox,N{:});
            obj.Type   = 'C2Q4';
        elseif strcmp(varargin{ii},'Hexahedron')
            N = num2cell(varargin{ii+1});
            obj.Center = Hexahedron(obj.BdBox,N{:});
            obj.Type   = 'C3H8';
        else
            obj.(varargin{ii}) = varargin{ii+1};
        end
    end
    
    if obj.Triangulate,  obj.Type = 'C2T3'; end
    if isempty(obj.SDF), obj.Type = 'C3T3'; end
    
    obj.Dim = 0.5*numel(obj.BdBox);
end

%---------------------------------------------------------------------- get     
function varargout = get(Mesh,varargin)
% gets any variable(s) in Class hierarchy.
% example usage: 
%  
%   [n,m] = msh.get('NNode','NElem')

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
% sets any variable(s) in Class hierarchy.
% example usage: 
%  
%   msh = msh.set('NNode',20,'NElem',42)
    for ii = 1:2:length(varargin)
        if strcmp(varargin{ii},'Quads')
            N = num2cell(varargin{ii+1});
            Mesh.Center = Quads(Mesh.BdBox,N{:});
            Mesh.Type = 'C2Q2';
        else
            Mesh.(varargin{ii}) = varargin{ii+1};
        end
    end
end

%---------------------------------------------------------------------- set
function Mesh = generate(Mesh)
    
Mesh.Iteration   = 1; 
Mesh.Convergence = 1e7;

if isempty(Mesh.Center)
    Pc = randomPointSet(Mesh); 
elseif strcmp(Mesh.Type,'C3H8')    
    Mesh.MaxIteration = 1;
    Pc = Mesh.Center; 
    Mesh.NElem = length(Pc);
else
    Mesh.MaxIteration = 1;
    Pc = Mesh.Center; 
    d = Mesh.SDF(Pc);
    Pc = Pc(d(:,end)<0,:);  
    Mesh.NElem = length(Pc);
end

% a-prioiri area estimate
if numel(Mesh.BdBox) < 6 
Anew = (Mesh.BdBox(2)-Mesh.BdBox(1))*(Mesh.BdBox(4)-Mesh.BdBox(3));
else
Anew = (Mesh.BdBox(2)-Mesh.BdBox(1))*(Mesh.BdBox(4)-Mesh.BdBox(3))*...
       (Mesh.BdBox(6)-Mesh.BdBox(5));
end

flag = 0;
% loyd's algorithm 
while flag == 0
  
  % update seeding points
  P = Pc;
  
  % reflection of seedings
  Rp = pointSetReflect(Mesh,P,Anew); 
  
  % bounding box seedings 
  Rb = boundingReflect(Mesh);

  % generate Voronoi tesselation
  [v,f] = voronoin([P;Rp;Rb],{'Qt','Pp'});   
  
  % compute centroid and area
  [Pc,A] = computeCentroid(Mesh,f,v); 
  
  % compute velocity
  Mesh.Velocity = vecnorm((P - Pc)')';
  
  Anew = sum(abs(A));
  
  Mesh.Convergence = vappend(Mesh.Convergence,sqrt(sum((A.^2).*sum((Pc-P)...
      .*(Pc-P),2)))*Mesh.NElem/(Anew^1.5));
  
  [flag,Mesh] = CheckConvergence(Mesh);
  
  Mesh.Iteration = Mesh.Iteration + 1;
  
  if Mesh.Movie
     Mesh.Center = Pc;
     Mesh.Node = v;
     Mesh.Element = f(1:Mesh.NElem);
     Mesh = show(Mesh,'Velocity');
  end
  
end

f = f(1:Mesh.NElem);

% %
[v,f] = RemoveDuplicates(Mesh,v,f);
[v,f] = ExtractNode(Mesh,v,f);
if Mesh.Dim < 3
[v,f] = CollapseEdges(Mesh,v,f); 
[v,f] = ResequenceNodes(Mesh,v,f);
end
if strcmp(Mesh.Type,'C3H8'), [v,f] = HexahedronOrder(Mesh,v,f); end

if Mesh.Triangulate
    Mesh.Node = v;
    Mesh.Element = f;
    [v,f] = MeshTriangulation(Mesh,P,v,f,length(f), length(v));
    [v,f] = ResequenceNodes(Mesh,v,f);
    Mesh.NElem = length(f);
    Mesh.NNode = length(v);
end

%[Pc,A,Nv,Ev] = computeCentroid(Mesh,f,v); 
[Pc,A] = computeCentroid(Mesh,f,v); 

Mesh.Center  = Pc;
Mesh.Node    = v;
Mesh.Element = f;
Mesh.NNode   = length(v);
Mesh.NElem   = length(f);
Mesh.Area    = A;
% Mesh.Normal  = Nv;
% Mesh.Edge  = Ev;

Mesh = ElementAdjecency(Mesh);
end

%---------------------------------------------------------------- show mesh
function Mesh = show(Mesh,varargin)
if nargin<2, Request = -1; 
else, Request = varargin{1}; end
figure(101);

% generate elemental matrices for plotting
Mesh = ElementAdjecency(Mesh);
fs = 'flat';

switch(Request)
    case('SDF'),      Z = Mesh.SDF(Mesh.Node); Z = Z(:,end);
    case('Velocity'), Z = (Mesh.Velocity')'; caxis([0 0.01]); 
    case('Gradient'), Z = GradientField(Mesh,Mesh.Node);
    case('Node'),     Z = varargin{2};
    case('Element'),  Z = varargin{2};
    otherwise;        Z = zeros(Mesh.NNode,1);
end

if ~strcmp(Request,'Node') && ~strcmp(Request,'Element')
clf; axis equal; axis off; hold on;
    
% plot tesselation
patch('Faces',Mesh.ElemMat,'Vertices',Mesh.Node,'FaceVertexCData',Z,...
    'Facecolor',fs,'LineStyle',Mesh.LineStyle,'Linewidth',1.5,...
    'EdgeColor','k');

% plot boundaries
patch('Faces',Mesh.Boundary,'Vertices',Mesh.Node,'LineStyle','-',...
    'Linewidth',2,'EdgeColor','k');

hold on;
elseif strcmp(Request,'Node')
    plot(Mesh.Node(Z,1),Mesh.Node(Z,2),'.','Color','r');
elseif strcmp(Request,'Element')
    plot(Mesh.Center(Z,1),Mesh.Center(Z,2),'.','Color','r');    
end

axis(Mesh.BdBox); axis tight;
colormap(Mesh.Colormap);
drawnow;

if Mesh.Movie 
    background(gitpage);
    if Mesh.MovieStart == false
       Mesh.MovieStart = true;
       MovieMaker(Mesh,'mesh','Start');
    else
       MovieMaker(Mesh,'mesh','');
    end
end

end

%--------------------------------------------------------------- class list
function list(Mesh)
% shows list of commands Mesh.list()
    
fprintf('* CLASS: Mesh() \n');
fprintf('\t :PUBLIC: \n');
fprintf('\t\t get(): \n');
fprintf('\t\t set(): \n');
fprintf('\t\t show(): \n');
fprintf('\t\t showSDF(): \n');
fprintf('\t\t generate(): \n');
fprintf('\t\t FindNodes(): \n');
fprintf('\t\t FindElements(): \n');

end

%----------------------------------------------------------------- show SDF
function showSDF(Mesh,varargin)
Bd = 1.05*Mesh.BdBox; 
Q  = 250;
xmin = Bd(1); xmax = Bd(2); 
ymin = Bd(3); ymax = Bd(4); 
%dl = 0.21*(0.5*((xmax - xmin) + (ymax - ymin))); 

% if xmin == 0, xmin = -Bd(2)/3; end
% if ymin == 0, ymin = -Bd(4)/3; end
x = linspace(xmin,xmax,Q);
y = linspace(ymin,ymax,Q);
[X,Y] = meshgrid(x,y);

P = [X(:),Y(:)];
D = Mesh.SDF(P);
Dist = reshape(D(:,end),[Q, Q]);
%Ds = reshape(mod(D(:,end),dl),[Q, Q]);
cla;
%image(x,y,Dist+2); hold on;
DistBnd = -wthresh(-Dist,'h',1e-6);


surf(x,y,Dist,'linestyle','none'); hold on;
contour3(X,Y,DistBnd,5,'w-','linewidth',1)%,%'Color',V(4,1:3)); 
contour(X,Y,Dist,[0 1e-6],'w-','linewidth',2); 
%colormap((bluesea(20)));
colormap(viridis);
caxis([-1,2 + 1e-6]);
axis equal;
axis off;
plotbox;
set(gca, 'YDir','normal');
view(0,90);
end

%--------------------------------------------------------------- find nodes
function NodeList = FindNodes(Mesh,Request)
    NodeList = FindNodes(Mesh.Node,Request);
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
P = zeros(Mesh.NElem,Mesh.Dim);
B = Mesh.BdBox;
Ctr = 0;
while(Ctr < Mesh.NElem)
    Y = zeros(Mesh.NElem,Mesh.Dim);
    for ii = 1:Mesh.Dim
    Y(:,ii) = (B(2*ii)-B(2*ii-1))*rand(Mesh.NElem,1)+B(2*ii-1);
    end
    d = Mesh.SDF(Y);
    I = find(d(:,end)<0);               
    NumAdded = min(Mesh.NElem-Ctr,length(I));
    P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
    Ctr = Ctr+NumAdded;
end
end 

%--------------------------------------------------------- reflect pointset
function Rp = pointSetReflect(Mesh,P,A)
Alpha = 1.5*(A/Mesh.NElem)^(1/2);
d = Mesh.SDF(P);  

% number of assigned boundary segments
NBdrySegs = size(d,2)-1;   
if Mesh.Dim == 2
n1 = (Mesh.SDF(P+repmat([Mesh.eps,0],Mesh.NElem,1))-d)/Mesh.eps;
n2 = (Mesh.SDF(P+repmat([0,Mesh.eps],Mesh.NElem,1))-d)/Mesh.eps;
else
n1 = (Mesh.SDF(P+repmat([Mesh.eps,0,0],Mesh.NElem,1))-d)/Mesh.eps;
n2 = (Mesh.SDF(P+repmat([0,Mesh.eps,0],Mesh.NElem,1))-d)/Mesh.eps;
n3 = (Mesh.SDF(P+repmat([0,0,Mesh.eps],Mesh.NElem,1))-d)/Mesh.eps;
end

% logical index of seeds near the boundary
I       = abs(d(:,1:NBdrySegs))<Alpha; 
P1      = repmat(P(:,1),1,NBdrySegs); 
P2      = repmat(P(:,2),1,NBdrySegs); 
Rp(:,1) = P1(I)-2*n1(I).*d(I);  
Rp(:,2) = P2(I)-2*n2(I).*d(I);
if Mesh.Dim == 3
P3 = repmat(P(:,3),1,NBdrySegs); 
Rp(:,3) = P3(I)-2*n3(I).*d(I);  
end

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
if Mesh.Dim == 3
tmp(5) = Mesh.BdBox(5)-a*( Mesh.BdBox(6) - Mesh.BdBox(5)); 
tmp(6) = Mesh.BdBox(6)+a*( Mesh.BdBox(6) - Mesh.BdBox(5)); 
end

if Mesh.Dim == 2
x = linspace(tmp(1),tmp(2),2);
y = linspace(tmp(3),tmp(4),2);
[X,Y] = meshgrid(x,y);
Rb = [X(:),Y(:)];  
else
x = linspace(tmp(1),tmp(2),2);
y = linspace(tmp(3),tmp(4),2);
z = linspace(tmp(5),tmp(6),2);
[X,Y,Z] = meshgrid(x,y,z);
Rb = [X(:),Y(:),Z(:)];      
end

end

%------------------------------------------------ compute centroid polygons
function [Pc,A] = computeCentroid(Mesh,f,v)
    
Pc = zeros(Mesh.NElem,Mesh.Dim); 
A  = zeros(Mesh.NElem,1);

if Mesh.Dim == 3
F = f(1:Mesh.NElem)';     
VCell = cellfun(@(E) v(E,:),F,'UniformOutput',false);
TCell  = cellfun(@(E) convhulln(E,{'Qt','QbB','Pp'}),VCell,...
'UniformOutput',false);

[Pc, A] = cellfun(@(V,F) CentroidPolyhedron(Mesh,V,F),VCell,...
TCell,'UniformOutput',false);

Pc = vertcat(Pc{:}); 
A = vertcat(A{:});

return;
end

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

%-------------------------------------- compute centroid convex poly-hedron
function [Pc,W] = CentroidPolyhedron(~,Node,Element)
    
    %finding the area of closed polyhedron
    v1 = Node(Element(:,2),:) - Node(Element(:,1),:);
    v2 = Node(Element(:,3),:) - Node(Element(:,1),:);
    
    Atmp = 0.5*cross(v1,v2);
    A(:,1) = sqrt(Atmp(:,1).^2+Atmp(:,2).^2+Atmp(:,3).^2);
    
    W = sum(A); %total area
    
    p1 = Node(Element(:,1),:);
    p2 = Node(Element(:,2),:);
    p3 = Node(Element(:,3),:);
    
    Pc = (1/3).*(p1 + p2 + p3); %centroid of each triangle
    
    mg(:,1) = A.*Pc(:,1);
    mg(:,2) = A.*Pc(:,2);
    mg(:,3) = A.*Pc(:,3);
    Mg = [mg(:,1),mg(:,2),mg(:,3)];
    
    Pc = sum(Mg)./W;
    
end

%----------------------------------------------- triangulate polygonal mesh
function [v, f] = MeshTriangulation(~,Center,v0,f0,Nel,Nvr)
f = [];

for ii = 1:Nel
    el = f0{ii};
    n = numel(el);
    elem = [el(1:n)', [el(2:n)'; el(1)], ...
            repmat(Nvr+ii,n,1)];
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

%------------------------------------------------------- extract nodal data
function [Node,Element] = RemoveDuplicates(Mesh,Node0,Element0)
[~,~,cNode] = unique(Node0,'rows');%uniquetol(Node0,'ByRows',1e-3);

[Node,Element] = Rebuild(Mesh,Node0,Element0(1:Mesh.NElem),cNode');
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

%---------------------------------------------------------- resequence mesh
function [Node,Element] = HexahedronOrder(Mesh,Node0,Element0)
Element = cell(Mesh.NElem,1);
Node = Node0;

for i = 1:Mesh.NElem
    v = Node0(Element0{i},:);
    s = C3H8Order(v);
    Element{i} = Element0{i}(s);
end

function s = C3H8Order(points)

% n = 8;
% k = n-1;
% s = 1:8;
% while (k > 0)
%     l = 0;
%     for j = 1:k
%     if (points(j,3) > points(j+1,3) ...
%     || (points(j,3) == points(j+1,3) && points(j,1) > points(j+1,1)) ...
%     || (points(j,3) == points(j+1,3) && points(j,1) == points(j+1,1)...
%         && points(j,2) < points(j+1,2)) )
%         for m = 1:3
%             tmp = points(j,m);
%             points(j,m)   = points(j+1,m);
%             points(j+1,m) = tmp;
%         end
%         tmp = s(j);
%         s(j) = s(j+1);
%         s(j+1) = tmp;
%         l = j;
%     end
%     end
%     k = l;
% end
% 
% s = [s(2),s(3),s(4),s(1),s(6),s(7),s(8),s(5)];
s = 1:8;
XD = [-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1];
for k = 1:8
    fvs = XD(:,k);

    if(fvs(1) == -1), [~,ix] = mink(points(:,1),4);
    else, [~,ix] = maxk(points(:,1),4);
    end

    if(fvs(2) == -1), [~,iy] = mink(points(:,2),4);
    else, [~,iy] = maxk(points(:,2),4);
    end

    if(fvs(3) == -1), [~,iz] = mink(points(:,3),4);
    else, [~,iz] = maxk(points(:,3),4);
    end
    
    tmp = intersect(intersect(ix,iy),iz);
    s(k) = tmp;
end

end

end

%------------------------------------------------ generate elemental matrix
function [ElemMat,RawConnect,id] = GenerateElementalMatrix(Mesh)
El = Mesh.Element(1:Mesh.NElem)';   
if Mesh.Dim == 3
nodeset = cellfun(@(X) Mesh.Node(X,:),Mesh.Element,'UniformOutput',false);
faceEl = cellfun(@(X) convhulln(X,{'Qt','Pp'}),nodeset,'UniformOutput',false);
faceDis = cellfun(@(V,F,E) CollectCoplanar(V,F,E),...
           nodeset,faceEl,Mesh.Element,'UniformOutput',false);
IdEl = cellfun(@(X,Y) Y*ones(size(X,1),1),faceDis,...
        num2cell(1:Mesh.NElem,1).','UniformOutput',false);       
face = num2cell(vertcat(faceDis{:}),2);
id = vertcat(IdEl{:});
El = face(:)'; 
R = cell(Mesh.NElem,1);
for ii = 1:length(id)
    R{id(ii)} = [R{id(ii)},face{id(ii)}]; 
end
Raw = R(:); 
Raw = reshape(Raw,[],1); MaxN = max(cellfun(@(E) size(E,2),Raw));        
PadWNaN = @(E) [E, NaN(size(E,1), MaxN- size(E,2))]; 
RawConnect = cellfun(PadWNaN,Raw,'UniformOutput',false);
RawConnect = vertcat(RawConnect{:});   
end

if Mesh.Dim < 3, RawConnect = []; id = []; end

El = reshape(El,[],1); MaxN = max(cellfun(@(E) size(E,2),El));        
PadWNaN = @(E) [E, NaN(size(E,1), MaxN- size(E,2))]; 
ElemMat = cellfun(PadWNaN,El,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});       
end

%----------------------------------------------- generate elemental adjency
function Mesh = ElementAdjecency(Mesh)
if Mesh.Dim == 2
    face = Mesh.Element;
elseif Mesh.Dim == 3
    nodeset = cellfun(@(X) Mesh.Node(X,:),Mesh.Element,'UniformOutput',false);
    if ~strcmp(Mesh.Type,'C3T3')
        faceEl = cellfun(@(X) convhulln(X,{'Qt','Pp'}),nodeset,...
            'UniformOutput',false);
    else
        faceEl = Mesh.Element;
    end
    face = num2cell(vertcat(faceEl{:}),2);
end

edges = cellfun(@numel,face);
n = Mesh.NElem;    
p = max(cellfun(@max,face));    
MaxNVer = max(cellfun(@numel,face));     
if Mesh.Dim < 3
    tri = GenerateElementalMatrix(Mesh);
    Mesh.ElemMat = tri;
else
    if ~strcmp(Mesh.Type,'C3T3')
        [A,tri,triID] = GenerateElementalMatrix(Mesh);
    else
        A = vertcat(Mesh.Element{:});
        tri = vertcat(Mesh.Element{:});
        triID = (1:Mesh.NElem).';
    end
    Mesh.ElemMat = A;
    MaxNVer = max(cellfun(@numel,num2cell(tri,2)));
    edges = cellfun(@numel,num2cell(tri,2));
    p = max(cellfun(@max,num2cell(tri,2)));
    Mesh.ElementToFace = sparse(1:length(A),triID,1);
end

set = 1:n;

%if Mesh.Dim == 2
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
%end
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

%------------------------------------------------------ isotropic reduction
function [flag,Mesh] = CheckConvergence(Mesh)
Criteria = (Mesh.Convergence(end) > Mesh.ConvNorm);

if (Criteria && (Mesh.Iteration < Mesh.MaxIteration)), flag = 0;
else
    if (Mesh.Iteration > Mesh.MaxIteration), flag = 2; 
    elseif (Mesh.Iteration == 1), flag = 0;        
    elseif strcmp(Mesh.Type,'C3H8'), flag = 2;        
    else, flag = 1;
    end
end
if Mesh.ShowProcess; ProcessMonitor(Mesh); end
end

%---------------------------------------------------------- material filter 
function ProcessMonitor(Mesh)
if Mesh.Iteration == 1
fprintf(' Iter   | Residual  | Max. Svm  | Time | dt      | p   |\n');
fprintf('--------------------------------------------------------------\n');
end
fprintf(' %i  \t| %1.3e | %1.3e |\n',...
    Mesh.Iteration,Mesh.Convergence(end),max(Mesh.Velocity));   
end
end

end

%-------------------------------------------------------------- movie maker
function MovieMaker(Mesh,Name,Request)
if nargin < 2, Request = ''; end

if Mesh.Movie
    switch(Request)
        case('Start')
            filename = string([Name,'_', char(datetime(now,...
                              'ConvertFrom','datenum')),'.gif']);
            
            filename = erase(filename,[":"," "]);
            background(gitpage); drawnow;
            gif(char(filename),'frame',gcf,'nodither','Timestep',1/12);
        otherwise
            background(gitpage);
            drawnow;
            gif;
    end
end

end

function R = PlanarProjection(a)
a = a(:); b = [0;0;1];
v = cross(a,b); vs = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
R = eye(3,3) + vs + vs*vs/(1+abs(dot(a,b))+1e-9);
end

%------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
function Poly = CollectCoplanar(Node,Element,Element0)
%finding the vector span of triangle
N = zeros(length(Element),3);
for i = 1:length(Element)
    v1 = Node(Element(i,2),:) - Node(Element(i,1),:);
    v2 = Node(Element(i,3),:) - Node(Element(i,1),:);   
    tmp = cross(v1,v2); N(i,:) = tmp/sqrt(tmp(1)^2 + tmp(2)^2 + tmp(3)^2);
end

[Norm,Ia] = uniquetol(N,1e-1,'ByRows',true,'OutputAllIndices', true);
Nc = transpose(num2cell(Norm.',1));

Rprj = cellfun(@(E) PlanarProjection(E),Nc,'UniformOutput',false); 
NodeId = cellfun(@(E) unique(Element(E,:)),Ia,'UniformOutput',false); 

Prj = cellfun(@(E,R)  R*transpose(Node(E,:)),NodeId,Rprj,'UniformOutput',false);
Stitch = @(E) convhull(transpose(E(1:2,:)));

Poly = cellfun(@(E) Stitch(E),Prj,'UniformOutput',false);
Poly = cellfun(@(E,V) reshape(E(V),1,[]), NodeId, Poly,'UniformOutput',false);
Poly = cellfun(@(E) Element0(E(:)), Poly,'UniformOutput',false);

Poly = NaNPadding(Poly);
Poly = vertcat(Poly{:}); 
end

%------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
function A = NaNPadding(A0)

A0 = reshape(A0,[],1); MaxN = max(cellfun(@(E) size(E,2),A0));        
PadWNaN = @(E) [E, NaN(size(E,1),MaxN- size(E,2))]; 
A = cellfun(PadWNaN,A0,'UniformOutput',false);

end

