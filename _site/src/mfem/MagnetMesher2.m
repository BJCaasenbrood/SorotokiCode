%  - structure - 
% ======= MagnetMesher.m =======
%  -- func.LocalElement
%  -- func.GenerateSetPoints
%  -- func.ElementDistribution
%  -- func.VoronoiLine
%  -- func.ComputeCentroid
% ===========================

function mesh = MagnetMesher(mesh,Domain)

mesh.BdBox = Domain('BdBox');
MaxIter = mesh.MaxIter;
Tol = 1e-6;       % tolerance

CallDisplay('mesh settings');
CallExecuted(['elements = ',num2str(mesh.NElem/1e3,2),'k']);
CallExecuted(['steps = ',num2str(mesh.MaxIter)]);

if mesh.Accel, CallExecuted('acceleration on');
else, CallExecuted('acceleration off'); end

if ~isfield(mesh,'P')
    [Pc, PFix, Sigma] = GenerateSetPoints(mesh,Domain);
else
    Pc = mesh.P;
    PFix = mesh.PFix;
end

% asking permission to mesh
if mesh.Permission, if ~Permission, return; end; end

CallDisplay('generating mesh...');

% while loop init
figure(1); Iter = 1; Err = Inf;

while(Iter <= MaxIter && Err > Tol)
    
    % Lloyd's update
    P = Pc;
    
    % make voronoi diagram
    [Pmid, Node, Element] = VoronoiLine(P,PFix,Domain);
    
    % Update centroid elements
    Pc = ComputeCentroid(Pmid,PFix,Domain);
    
    if mesh.Accel ~= 1
        mesh = WriteMesh(mesh,Node,Element);
        PlotMagnetMesh(mesh);
        fprintf('It: %3d   Error: %1.3e\n',Iter,Err);
    else
        %mesh = WriteMesh(mesh,Node,Element);
    end
    
    Err = sqrt(sum((Pc-P).^2)); 
    Error(Iter) = sqrt(sum(Pc-P).^2);
    Iter = Iter + 1;
end

Node = Domain('Curve',sort([Pc(:);PFix(:)]));
mesh = WriteMesh(mesh,Node,Element);

mesh.FigHandle = PlotMagnetMesh(mesh);
mesh.Error = Error;
mesh.Current = Domain('BC');

CallExecuted('meshing done!')

end

%----------------------------------------------- generate spatial setpoints
function [P,PFix,Sigma]  = GenerateSetPoints(mesh,Domain)
Sgm = Domain('Segments');
NElem = mesh.NElem;

% generate 1D-spatial domains
D = [0,1]; E = zeros(Sgm,2);
for i = 2:Sgm, D(i,:) =  D(1,:) - D(1,:)/2^(i-1); end; D = unique(D);
for i = 1:Sgm, E(i,1) = D(i); E(i,2) = D(i+1); end

Ctr = ElementDistribution(Domain,NElem,E);

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
function Ctr = ElementDistribution(Domain,NElem,E)

Sgm = Domain('Segments'); 
N = 1e5;
Vol = zeros(Sgm,1);

for i = 1:Sgm
    P = (E(i,2)-E(i,1))*transpose(linspace(1e-6,1-1e-6,N))+E(i,1);
    C = Domain('Curve',P(:)); 
    dC = num2cell(diff(C),2);
    dL = cellfun(@(E) norm(E,2),dC,'UniformOutput',false);
    Vol(i) = sum(vertcat(dL{:}));
end

L = Vol/sum(Vol); Ctr = zeros(Sgm,1);
for i = 1:Sgm-1, Ctr(i) = ceil(NElem*L(i)); end

Ctr(end) = NElem - sum(Ctr);

end

%------------------------------------------------ generate linear setpoints
function [Pmid, Node, Element] = VoronoiLine(P,PFix,Domain)

for j = 1:length(PFix)-1
    id = P>PFix(j) & P<PFix(j+1);
    SubP{j} = P(id);%[PFix(j);P(id);PFix(j+1)];
    subMidP{j} = diff(SubP{j})/2 + SubP{j}(1:end-1);
end

Pmid = vertcat(subMidP{:});
PNew = sort([Pmid;PFix(:)]);

Node = Domain('Curve',PNew(:));

E0 = 1:length(PNew); E1 = E0+1;
Element = [E0(:),E1(:)];
Element(end,:) = [];
end

%------------------------------------------------ generate linear setpoints
function Pc = ComputeCentroid(P,PFix,Domain)

for j = 1:length(PFix)-1
    id = P>=PFix(j) & P<=PFix(j+1);
    subP{j} = [PFix(j);P(id);PFix(j+1)];
    subPCell{j} = [subP{j}(1:end-1),subP{j}(2:end)];
end

PCells = num2cell(vertcat(subPCell{:}),2);
     
dX = cellfun(@(E) VirtualCenter(E,Domain),PCells,'UniformOutput',false);
         
%Pc = P + vertcat(dX{:});
Pc = vertcat(dX{:});

end

%---------------------------------------------------------------- PLOT MESH
function P = VirtualCenter(Pfield,Domain)
N = 3;
dp = linspace(Pfield(1),Pfield(2),N);
dLc = diff(Domain('Curve',dp(:)));
dL = cumsum(sqrt(dLc(:,1).^2 + dLc(:,2).^2 + dLc(:,3).^2));

[~,id] = min(abs(dL - dL(end)/2));

P = dp(id+1);
end

%---------------------------------------------------------------- PLOT MESH
function I = ComputeCurrentVectors(Node,Element)

NNode = length(Node);
I = zeros(NNode,3);
for ii = 1:NNode-1
    N1 = Node(ii,:);
    N2 = Node(ii+1,:);
    I(ii,:) = (N1 - N2)/norm(N1 - N2);
    I(ii,1) = 0.5*I(ii,1) + 0.5;
    I(ii,2) = 0.5*I(ii,2) + 0.5;
    I(ii,3) = 0.5*abs(I(ii,3)) + 0.5;
end

I(end,:) = I(end-1,:);

end

%---------------------------------------------------------------- PLOT MESH
function FigHandle = PlotMagnetMesh(mesh)

BdBox = mesh.BdBox;
Node = mesh.Node;
Element = mesh.Element;

clf; axis equal; axis off; hold on;

FigHandle = patch('Faces',Element,'Vertices',Node,'LineStyle','-',...
    'Edgecolor',[233, 84, 32]/255,'Linewidth',1,...
    'Marker','none','Markersize',5);

PlotGroundPlane(BdBox);
PlotAxes(BdBox);

axis(1.1*BdBox); view(110,20); axis equal; axis vis3d;
set(gcf,'color',[255, 255, 255]/255);
pause(1e-6);
axis on;
drawnow;

end


%------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
function PlotAxes(BdBox)

tmp = BdBox; 

a = 0.2;
tmp(1) = BdBox(1)-a*( BdBox(2) - BdBox(1)); 
tmp(2) = BdBox(2)+a*( BdBox(2) - BdBox(1)); 
tmp(3) = BdBox(3)-a*( BdBox(4) - BdBox(3)); 
tmp(4) = BdBox(4)+a*( BdBox(4) - BdBox(3)); 
tmp(5) = BdBox(5)-a*( BdBox(6) - BdBox(5)); 
tmp(6) = BdBox(6)+a*( BdBox(6) - BdBox(5)); 
BdBox = tmp;

X = [BdBox(1),BdBox(3),0;BdBox(2),BdBox(3),0]; dX = diff(X)*(1+0.52*a);
Y = [BdBox(1),BdBox(3),0;BdBox(1),BdBox(4),0]; dY = diff(Y)*(1+0.52*a);
%Z = [BdBox(1),BdBox(3),0;BdBox(1),BdBox(3),BdBox(6)]; dZ = diff(Z)*(1+0.52*a);

quiver3(X(1,1),X(1,2),X(1,3),dX(1,1),dX(1,2),dX(1,3),'r','Linewidth',2);
quiver3(Y(1,1),Y(1,2),Y(1,3),dY(1,1),dY(1,2),dY(1,3),'g','Linewidth',2);
%quiver3(Z(1,1),Z(1,2),Z(1,3),dZ(1,1),dZ(1,2),dZ(1,3),'b','Linewidth',2);

end

%--------------------------------------------------------- PLOT GROUNDPLANE
function PlotGroundPlane(BdBox)
tmp = BdBox; a = 0.2;
tmp(1) = BdBox(1)-a*( BdBox(2) - BdBox(1));
tmp(2) = BdBox(2)+a*( BdBox(2) - BdBox(1));
tmp(3) = BdBox(3)-a*( BdBox(4) - BdBox(3));
tmp(4) = BdBox(4)+a*( BdBox(4) - BdBox(3));
tmp(5) = BdBox(5)-0*( BdBox(6) - BdBox(5));
tmp(6) = BdBox(6)+a*( BdBox(6) - BdBox(5));
BdBox = tmp;

Nx = 4;
Ny = 4;
DN = 2;
x = linspace(BdBox(1),BdBox(2),Nx+1);
y = linspace(BdBox(3),BdBox(4),Ny+1);
xh = linspace(BdBox(1),BdBox(2),(DN+1)*(Nx)-(Nx-1));
yh = linspace(BdBox(3),BdBox(4),(DN+1)*(Ny)-(Ny-1));

[X,Y] = meshgrid(x,y);
[Xh,Yh] = meshgrid(xh,yh);

offset = -1e-5*(BdBox(6)-BdBox(5)) + BdBox(5);

mesh(Xh,Yh,Xh*0+0*offset,'Edgecolor',ColorScheme(1),...
     'Facecolor','none','linewidth',.1);

mesh(X,Y,X*0+0*offset,'Edgecolor',ColorScheme(1),...
     'Facecolor','none','linewidth',1.5);

end

%------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
function mesh = WriteMesh(mesh,Node,Element)

mesh.Node = Node;
mesh.Element = Element;
mesh.NElem = size(Element,1);
mesh.NNode = size(Node,1);

assignin('base','mesh',mesh);

end

%------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
function bool = Permission

request = CallRequest('initiate meshing?','y/n');

if strcmp(request,'n'), CallWarning('meshing terminated','Mesher3D:'); bool = false;
elseif strcmp(request,'y'), bool = true;
else, CallWarning('meshing terminated','Mesher3D:');  bool = false;
end

end
