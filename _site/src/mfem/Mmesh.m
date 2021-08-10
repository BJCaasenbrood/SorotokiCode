classdef Mmesh

    properties (Access = public)
        SDF;
        BdBox;
        Node;
        Node0;
        Element;
        Segment;
        NSegment;
        NNode;
        NElem;
        Dim;
        Center;
        Anchors;
        Inductance;
    end
    
    properties (Access = private)
        Current;
        Countbin;
        Area;
        WireThickness;
        Convergence;
        Velocity;
        Adjecency;
        Normal;
        Edge;
        Iteration;
        ElemMat;
        eps; 
        eta; 
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
    
function obj = Mmesh(Input,varargin) 
    
    if isa(Input,'function_handle')
       obj.SDF = Input;
       obj.NElem = 200;
    end

    obj.NSegment     = 1;
    obj.MaxIteration = 100;
    obj.Iteration    = 0;
    obj.ElemMat      = -1;
    obj.eps          = 1e-6; 
    obj.eta          = 0.9; %was 0.9
    obj.Movie        = false;
    obj.MovieStart   = false;
    obj.CollapseTol  = 0.2;
    obj.ConvNorm     = 1e-3;
    obj.ShowProcess  = false;
    obj.Colormap     = turbo;
    obj.LineStyle    = '-';
    obj.WireThickness = 0.05;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
end

%---------------------------------------------------------------------- get     
function varargout = get(Mmesh,varargin)
% gets any variable(s) in Class hierarchy.
% example usage: 
%  
%   [n,m] = msh.get('NNode','NElem')
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Mmesh.(varargin{ii});
        end
    else
        varargout = Mmesh.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function Mmesh = set(Mmesh,varargin)
% sets any variable(s) in Class hierarchy.
% example usage: 
%  
%   msh = msh.set('NNode',20,'NElem',42)
    for ii = 1:2:length(varargin)
       Mmesh.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- set
function Mmesh = generate(Mmesh)
    
Mmesh.Iteration   = 1; 
Mmesh.Convergence = 1e7;

if isempty(Mmesh.Center)
    [Pc, PFix, Sigma] = RandomPointSet(Mmesh);
    Mmesh.Anchors  = PFix;
    Mmesh.Countbin = Sigma;
else
    Pc   = Mmesh.Center;
    %PFix = Mmesh.Anchors;
end

flag = 0;
% loyd's algorithm 
while flag == 0
  
  % update seeding points
  P = Pc;

  % generate Voronoi tesselation
  [p, v, f] = VoronoiLine(Mmesh,[P;PFix]);
  
  % compute centroid and area
  Pc = ComputeCentroid(Mmesh,p);
  
  % compute velocity
  Mmesh.Velocity = vecnorm((P - Pc)')';
  
  % compute convergence
  Mmesh.Convergence = vappend(Mmesh.Convergence,sqrt(sum(sum((Pc-P)...
      .*(Pc-P),2)))*Mmesh.NElem);
  
  [flag,Mmesh] = CheckConvergence(Mmesh);
  
  Mmesh.Iteration = Mmesh.Iteration + 1;
  
end

v = Mmesh.SDF(sort([Pc(:);PFix(:)]));

Mmesh.Node    = v;
Mmesh.Node0   = v;
Mmesh.Element = f;
Mmesh.NNode = length(v);
Mmesh.NElem = length(f);
Mmesh.BdBox = boxhull(v);

if abs( Mmesh.BdBox(2) - Mmesh.BdBox(1)) < 1e-2
    Mmesh.BdBox(1) = Mmesh.BdBox(3);
    Mmesh.BdBox(2) = Mmesh.BdBox(4);
elseif abs( Mmesh.BdBox(4) - Mmesh.BdBox(3)) < 1e-2
    Mmesh.BdBox(3) = Mmesh.BdBox(1);
    Mmesh.BdBox(4) = Mmesh.BdBox(2);
end

C = [1;cumsum(Mmesh.Countbin)];

% C = Mmesh.SDF(PFix);
% [~,NFix] = ismember(C,v,'rows');
% 
for ii = 1:length(C)-1
    Mmesh.Segment{ii,1} = ((C(ii)):(C(ii+1))).';
end

Mmesh = ComputeCurrentVectors(Mmesh);

end

%---------------------------------------------------------------------- set
function [Mmesh,Y] = inductance(Mmesh,varargin)
   
Mmesh.Iteration = 1; 
Rmin = 1e-6;
mu0 = 4*pi*1e-7;
dcf = Mmesh.NElem*1e3;
Ltmp = zeros(Mmesh.NElem,Mmesh.NElem); 
    
if isempty(varargin)
    [dC,dR,Cm] = LocalElement(Mmesh.Node,Mmesh.Element);
    PermMat = PermutationMatrix(Mmesh.NElem);
else
    id = Mmesh.Segment{varargin{1}};
    [dC,dR,Cm] = LocalElement(Mmesh.Node,Mmesh.Element(id,:));
    PermMat = PermutationMatrix(length(Mmesh.Element(id,:)));
end

if ~PermMat, return; end
    
for k = 1:length(PermMat)
    id = PermMat(k,1);
    jd = PermMat(k,2);
    
    Rtmp = norm(Cm(id,:)-Cm(jd,:));
    
    % check if 
    if Rtmp > Rmin/2
        Ltmp(id,jd) = dC(id,:)*dC(jd,:).'/Rtmp;
        Ltmp(jd,id) = Ltmp(id,jd);
    end
    
    if ~mod(Mmesh.Iteration,dcf)
        fprintf(['  -batch: ',num2str((Mmesh.Iteration-dcf)/1e3),'k','-',...
            num2str(Mmesh.Iteration/1e3),'k\n']);
    end
    
    Mmesh.Iteration = Mmesh.Iteration + 1;
end

Y = mu0*sum(sum(Ltmp))/(4*pi);
Mmesh.Inductance = Y;
cprintf('green','* inductance = %1.3e \n',Y);

    
end

%--------------------------------------------------------------------- show
function FigHandle = show(Mmesh)

B = Mmesh.BdBox;

clf; axis equal; axis off; hold on;

% [x,y,z] = tubeplot(Mmesh.Node.',3e-4);
% h = surf(x,y,z);
% fv = surf2patch(h);

% FigHandle = patch('Faces',fv.faces,'Vertices',fv.vertices,'LineStyle','-',...
%       'FaceVertexCData',Mmesh.Current,'Linewidth',1,'Edgecolor','interp',...
%     'Marker','none','Markersize',5);

FigHandle = patch('Faces',Mmesh.Element,'Vertices',Mmesh.Node,'LineStyle','-',...
      'FaceVertexCData',Mmesh.Current,'Linewidth',1,'Edgecolor','interp',...
    'Marker','none','Markersize',5);

% PlotGroundPlane(BdBox);
% PlotAxes(BdBox);
% 
% axis(1.1*BdBox); view(110,20); axis equal; axis vis3d;
% set(gcf,'color',[255, 255, 255]/255);
colormap(Mmesh.Colormap);
% pause(1e-6);
% axis on;
drawnow;

end

%--------------------------------------------------------------------- show
function render(Mmesh)
    
    class = whoClasses('Gmodel');
    for i = 1:length(class)
        if ~isempty(class{i})
            delete(class{i}.get('FigHandle'));
        end
    end
    
    obj = Gmodel(Mmesh);
    obj.Texture = copper;
    
    obj.bake().render();

    assignin('base','obj',obj);
    
end

%-------------------------------------------------------------- plot ground
function Mmesh = ground(Mmesh,gnd)
    if nargin < 2, Groundplane(Mmesh);
    else, Groundplane(Mmesh,gnd); end
end

%-------------------------------------------------------------- plot ground
function Mmesh = box(Mmesh)
    BoundingBox(Mmesh);
end

%-------------------------------------------------------------- plot ground
function Mmesh = reset(Mmesh)
Mmesh.Node = Mmesh.Node0;
end

%-------------------------------------------------------------- END METHODS
end
methods (Access = private)

%----------------------------------------------- generate spatial setpoints
function [P,PFix,Sigma]  = RandomPointSet(Mmesh)
Sgm = Mmesh.NSegment;
NEl = Mmesh.NElem-Sgm;

% generate 1D-spatial domains
D = [0,1]; E = zeros(Sgm,2);
for i = 2:Sgm, D(i,:) =  D(1,:) - D(1,:)/2^(i-1); end; D = unique(D);
for i = 1:Sgm, E(i,1) = D(i); E(i,2) = D(i+1); end

Ctr = ElementDistribution(Mmesh,NEl,E);

Sigma = 0;
while(min(Sigma)==0)
    
P = [];
for i = 1:Sgm
    P = [P;(E(i,2)-E(i,1))*rand(Ctr(i),1)+E(i,1)];
end

P = sort(P(:)); 
PFix = D(:); 
Sigma = Ctr;

end

end

%------------------------------------------------ generate linear setpoints
function Ctr = ElementDistribution(Mmesh,NEl,E)

Sgm = Mmesh.NSegment;
N = 1e5;
Vol = zeros(Sgm,1);

for i = 1:Sgm
    P = (E(i,2)-E(i,1))*transpose(linspace(1e-6,1-1e-6,N))+E(i,1);
    C = Mmesh.SDF(P(:));
    dC = num2cell(diff(C),2);
    dL = cellfun(@(E) norm(E,2),dC,'UniformOutput',false);
    Vol(i) = sum(vertcat(dL{:}));
end

L = Vol/sum(Vol); Ctr = zeros(Sgm,1);
for i = 1:Sgm-1, Ctr(i) = ceil(NEl*L(i)); end

Ctr(end) = NEl - sum(Ctr);

end

%------------------------------------------------ generate linear setpoints
function [p, v, f] = VoronoiLine(Mmesh,P)
%Node = Domain('Curve',sort(P(:)));
PFix = Mmesh.Anchors;

for j = 1:length(PFix)-1
    id = P>PFix(j) & P<PFix(j+1);
    SubP{j} = P(id);%[PFix(j);P(id);PFix(j+1)];
    subMidP{j} = diff(SubP{j})/2 + SubP{j}(1:end-1);
end

p = vertcat(subMidP{:});
PNew = sort([p;PFix(:)]);

v = Mmesh.SDF(PNew(:)); 

E0 = 1:length(PNew); E1 = E0+1;
f = [E0(:),E1(:)];
f(end,:) = [];
end

%------------------------------------------------ generate linear setpoints
function Pc = ComputeCentroid(Mmesh,P)

PFix = Mmesh.Anchors;

for j = 1:length(PFix)-1
    id = P>=PFix(j) & P<=PFix(j+1);
    subP{j} = [PFix(j);P(id);PFix(j+1)];
    subPCell{j} = [subP{j}(1:end-1),subP{j}(2:end)];
end
% 
P = subPCell{:};
Pp = [P(:,1), (P(:,1) + P(:,2))/2, P(:,2)];

%KCells = mat2cell((Mmesh.SDF(Pp(:))),3*ones(length(P),1));
PCells = num2cell(vertcat(subPCell{:}),2); 
dX = cellfun(@(E) VC(E,Mmesh),PCells,'UniformOutput',false);
%dX = cellfun(@(E,P) VirtualCenter(E,P),KCells,PCells,'UniformOutput',false);
Pc = vertcat(dX{:});

% subfunction routine
function dx = VirtualCenter(K,P)
    dp = [P(1),(P(1) + P(2))/2,P(2)];
    dLc = diff(K);
    dL = cumsum(sqrt(dLc(:,1).^2 + dLc(:,2).^2 + dLc(:,3).^2));
 
    [~,k] = min(abs(dL - dL(end)/2));
    
    dx = dp(k+1);
end

function P = VC(Pfield,Mmesh)
N = 3;
dp = linspace(Pfield(1),Pfield(2),N);
dLc = diff(Mmesh.SDF(dp(:)));
dL = cumsum(sqrt(dLc(:,1).^2 + dLc(:,2).^2 + dLc(:,3).^2));

[~,k] = min(abs(dL - dL(end)/2));

P = dp(k+1);
end

end

%---------------------------------------------------------------- PLOT MESH
function Mmesh = ComputeCurrentVectors(Mmesh)

v = Mmesh.Node;
I = zeros(Mmesh.NNode,3);

for ii = 1:Mmesh.NNode-1
    N1 = v(ii,:);
    N2 = v(ii+1,:);
    I(ii,:) = (N1 - N2)/norm(N1 - N2);
    I(ii,1) = 0.5*I(ii,1) + 0.5;
    I(ii,2) = 0.5*I(ii,2) + 0.5;
    I(ii,3) = 0.5*abs(I(ii,3)) + 0.5;
end

%I(end,:) = I(end-1,:);
Mmesh.Current = I;
end

%---------------------------------------------------------------- PLOT MESH
function Element = RemoveStringElements(Mmesh)
Node = Mmesh.Node;
El = num2cell(Element,2);
Euclead = @(E) ((Node(E(1),1) - Node(E(2),1)).^2 + (Node(E(1),2) - ...
    Node(E(2),2)).^2 + (Node(E(1),3) - Node(E(2),3)).^2)^0.5;
dx = cellfun(@(E) Euclead(E),El);

end

%------------------------------------------------------ isotropic reduction
function [flag,Mmesh] = CheckConvergence(Mmesh)
Criteria = (Mmesh.Convergence(end) > Mmesh.ConvNorm);

if (Criteria && (Mmesh.Iteration < Mmesh.MaxIteration)), flag = 0;
else
    if (Mmesh.Iteration > Mmesh.MaxIteration), flag = 2; 
    elseif (Mmesh.Iteration == 1), flag = 0;           
    else, flag = 1;
    end
end

if Mmesh.ShowProcess; ProcessMonitor(Mmesh); end
end

%---------------------------------------------------------- material filter 
function ProcessMonitor(Mmesh)
if Mmesh.Iteration == 1
fprintf(' Iter   | Residual  | Max. vel  |\n');
fprintf('--------------------------------------------------------------\n');
end
fprintf(' %i  \t| %1.3e | %1.3e |\n',...
    Mmesh.Iteration,Mmesh.Convergence(end),max(Mmesh.Velocity));   
end

%----------------------------------------------------- generate groundplane
function Groundplane(Mmesh,gnd)
if nargin < 2
    tmp = Mmesh.BdBox; a = 0.1;    
    tmp(1) = Mmesh.BdBox(1)-a*( Mmesh.BdBox(2) - Mmesh.BdBox(1)); 
    tmp(2) = Mmesh.BdBox(2)+a*( Mmesh.BdBox(2) - Mmesh.BdBox(1)); 
    tmp(3) = Mmesh.BdBox(3)-a*( Mmesh.BdBox(4) - Mmesh.BdBox(3)); 
    tmp(4) = Mmesh.BdBox(4)+a*( Mmesh.BdBox(4) - Mmesh.BdBox(3)); 
    tmp(5) = Mmesh.BdBox(5)-0*( Mmesh.BdBox(6) - Mmesh.BdBox(5))-1e-3; 
    tmp(6) = Mmesh.BdBox(6)+a*( Mmesh.BdBox(6) - Mmesh.BdBox(5)); 
else
    tmp = gnd;
end

B = tmp;

Nx = 4;
Ny = 4;

x = linspace(B(1),B(2),Nx+1);
y = linspace(B(3),B(4),Ny+1);

[X,Y] = meshgrid(x,y);

v = [tmp(1),tmp(3), tmp(5);  
     tmp(2),tmp(3), tmp(5);
     tmp(2),tmp(4), tmp(5);
     tmp(1),tmp(4), tmp(5)];
 
f = [1,2,3,4];

[I,~] = imread('checker.jpg');

dX = (tmp(2)-tmp(1));
dY = (tmp(4)-tmp(3));

if dY/dX >= 2, I = vertcat(I,I);
elseif dX/dY >= 2, I = horzcat(I,I);
end

hold all
warpim(X,Y,X*0 + tmp(5),I);
hold off;

patch('Faces',f,'Vertices',v,...
    'Linewidth',1.5,'linestyle','-','FaceColor','none',...
    'EdgeColor',[1 1 1]*0.5);
end

%---------------------------------------------------- generate bounding box
function BoundingBox(Mmesh)

tmp = boxhull(Mmesh.Node);

v = [tmp(1),tmp(3), tmp(5);  
     tmp(2),tmp(3), tmp(5);
     tmp(2),tmp(4), tmp(5);
     tmp(1),tmp(4), tmp(5);
     tmp(1),tmp(3), tmp(6);  
     tmp(2),tmp(3), tmp(6);
     tmp(2),tmp(4), tmp(6);
     tmp(1),tmp(4), tmp(6)];
 
f = [1,2,3,4;
     5,6,7,8;
     1,5,nan,nan;
     2,6,nan,nan;
     3,7,nan,nan;
     4,8,nan,nan];

patch('Faces',f,'Vertices',v,...
    'Linewidth',1.5,'linestyle','-','FaceColor','none',...
    'EdgeColor',col(2));

end

end
end


% %---------------------------------------------------------------- PLOT MESH
% function FigHandle = PlotMagnetMesh(mesh)
% 
% BdBox = mesh.BdBox;
% Node = mesh.Node;
% NNode = mesh.NNode;
% Element = mesh.Element;
% Curr = mesh.Current;
% 
% clf; axis equal; axis off; hold on;
% 
% ColorCode = parula(NNode);
% 
% %ColorCode = @(x) repmat(NN,NNode,1) + 0.5*(x-1)*ColorScheme(1) + 0.5*x*ColorScheme(2);
% Color = [233, 84, 32]/255;
% 
% FigHandle = patch('Faces',Element,'Vertices',Node,'LineStyle','-',...
%     'FaceVertexCData',Curr,'Linewidth',1,'Edgecolor','interp',...
%     'Marker','none','Markersize',5);
% 
% PlotGroundPlane(BdBox);
% 
% axis(1.1*BdBox); view(110,20); axis equal; axis vis3d;
% set(gcf,'color',[255, 255, 255]/255);
% pause(1e-6);
% axis off;
% drawnow;
% 
% end
% 
% %------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
% function mesh = WriteMesh(mesh,Node,Element,Pc)
% 
% mesh.Node = Node;
% mesh.Element = Element;
% mesh.NElem = size(Element,1);
% mesh.NNode = size(Node,1);
% mesh.Center = Pc;
% 
% assignin('base','mesh',mesh);
% 
% end
% 
% %------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
% function bool = Permission
% 
% request = CallRequest('initiate meshing?','y/n');
% 
% if strcmp(request,'n'), CallWarning('meshing terminated','Mesher3D:'); bool = false;
% elseif strcmp(request,'y'), bool = true;
% else, CallWarning('meshing terminated','Mesher3D:');  bool = false;
% end
% 
% end
% end
