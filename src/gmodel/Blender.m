function x = Blender(x,Request,Arg)

switch(Request)
    case('LoadMesh');     x = LoadMesh(Arg);
    case('Translate');    x = TranslateMesh(x,Arg);
    case('Rotate');       x = RotateMesh(x,Arg);
    case('Center');       x = CenterMesh(x);
    case('Transform');    x = Transformation(x,Arg);
    case('Scale');        x = ScaleMesh(x,Arg);
    case('Curve');        x = CurveMesh(x,Arg);
    case('Sweep');        x = SweepMesh(x,Arg);
    case('SE3');          x = SE3Mesh(x,Arg);
    case('Twist');        x = TwistMesh(x,Arg);
end

end

%---------------------------------------------------------------- LOAD MESH
function mesh = LoadMesh(Arg)
mesh = struct;
path = Arg;
type = path(end-2:end);

if strcmp(type,'stl')
    [Node, Element] = stlreader(path);
elseif strcmp(type,'obj')
    [Node, Element] = objreader(path);
else
    error('File type not recognized or supported!');
end

mesh.Node = Node;
mesh.Element = Element;
mesh.BdBox = [min(Node(:,1)), max(Node(:,1)),...
              min(Node(:,2)), max(Node(:,2)),...
              min(Node(:,3)), max(Node(:,3))];
end

%--------------------------------------------------------------- SCALE MESH
function mesh = TranslateMesh(mesh,Arg)

Ax = Arg{1};
Move  = Arg{2};
Node0 = mesh.Node;
Node = Node0;

if strcmp(Ax,'x'),     Node(:,1) = Node0(:,1) + Move;
elseif strcmp(Ax,'y'), Node(:,2) = Node0(:,2) + Move;
elseif strcmp(Ax,'z'), Node(:,3) = Node0(:,3) + Move;
elseif strcmp(Ax,'3D')
    dx  = Arg{2}; Node(:,1) = Node0(:,1) + dx;
    dy  = Arg{3}; Node(:,2) = Node0(:,2) + dy;
    dz  = Arg{4}; Node(:,3) = Node0(:,3) + dz;  
else, Node(:,3) = Node0(:,3) + Move;
end

mesh.Node = Node;

end

%-------------------------------------------------------------- ROTATE MESH
function mesh = RotateMesh(mesh,Arg)
Ax = Arg{1};
dr  = Arg{2}*(pi/180);

Node0 = mesh.Node;

if strcmp(Ax,'x')
    R = [1,0,0; 0, cos(dr),-sin(dr); 0, sin(dr), cos(dr)];
    Node = transpose(R*(Node0'));
elseif strcmp(Ax,'y')
    R = [cos(dr), 0,-sin(dr);0,1,0;sin(dr), 0, cos(dr)];
    Node = transpose(sparse(R)*(Node0'));
elseif strcmp(Ax,'z')
    R = [cos(dr),-sin(dr), 0; sin(dr), cos(dr), 0; 0,0,1];
    Node = transpose(sparse(R)*(Node0'));
elseif strcmp(Ax,'3D')
    drx  = Arg{2}*(pi/180);
    dry  = Arg{3}*(pi/180);
    drz  = Arg{4}*(pi/180);
    c_3 = cos(drx); s_3 = sin(drx);
    c_2 = cos(dry); s_2 = sin(dry);
    c_1 = cos(drz); s_1 = sin(drz);
    rotm = zeros(3,3);
    rotm(1,1) =  c_1*c_2;
    rotm(1,2) =  c_1*s_2*s_3 - s_1*c_3;
    rotm(1,3) =  c_1*s_2*c_3 + s_1*s_3;
    
    rotm(2,1) =  s_1*c_2;
    rotm(2,2) =  s_1*s_2*s_3 + c_1*c_3;
    rotm(2,3) =  s_1*s_2*c_3 - c_1*s_3;
    
    rotm(3,1) = -s_2;
    rotm(3,2) =  c_2*s_3;
    rotm(3,3) =  c_2*c_3;
    
    Node = transpose(rotm*(Node0'));
end

mesh.Node = Node;

end

%-------------------------------------------------------------- ROTATE MESH
function Node = Transformation(mesh,Arg)

Node0 = mesh.Node; 
Node0(:,4) = 1;

HMat = Arg;

Node = HMat*Node0.';
Node = Node(1:3,:).';
mesh.Node = Node;
end

%-------------------------------------------------------------- ROTATE MESH
function mesh = CenterMesh(mesh)

Node0 = mesh.Node; 

X0 = mean(Node0(:,1));
Y0 = mean(Node0(:,2));
Z0 = mean(Node0(:,3));

mesh.Node = Node0 - [X0,Y0,Z0];
end

%-------------------------------------------------------------- ROTATE MESH
function mesh = SE3Mesh(mesh,Arg)

Node0 = [mesh.Node, ones(mesh.NNode,1)]; 

R = Quat2Rot(Arg(1:4));
r = Arg(5:7);
H = zeros(4);
H(1:3,1:3) = R.';
H(1:3,4) = [r(3),r(2),r(1)].';
H(4,4) = 1;

Node = H*Node0.';
mesh.Node = Node(1:3,:).';

end

%--------------------------------------------------------------- SCALE MESH
function mesh = ScaleMesh(mesh,Arg)
if isa(Arg,'cell')
    Ax = Arg{1};
    Scale  = Arg{2};
else
   Ax = 'uniform';
   Scale = Arg;
end
    
Node0 = mesh.Node;
Node = Node0;

if strcmp(Ax,'x')
    Node(:,1) = Scale*Node0(:,1);
    Node(:,1) = Node(:,1);% - mean(Node0(:,1));
elseif strcmp(Ax,'y')
    Node(:,2) = Scale*Node0(:,2);
    Node(:,2) = Node(:,2);% - mean(Node0(:,2));
elseif strcmp(Ax,'z')
    Node(:,3) = Scale*Node0(:,3);
    Node(:,3) = Node(:,3);% - min(Node0(:,3));
elseif strcmp(Ax,'axi')
    Node(:,1) = Scale*Node0(:,1);
    Node(:,1) = Node(:,1);  
    Node(:,2) = Scale*Node0(:,2);
    Node(:,2) = Node(:,2);
elseif strcmp(Ax,'uniform')
    Node = Scale*Node;
else 
end

mesh.Node = Node;

end

%--------------------------------------------------------------- SCALE MESH
function mesh = SweepMesh(mesh,Arg)

Node0 = mesh.Node;
List = Arg{1};
Swp = Arg{2};
V = 0*Node0;

for ii = 1:length(Node0)
    v = Node0(ii,:);
    vx = v(1);
    vy = v(2);
    id = List(ii);
    
    p1 = [vx; vy; 0];
    p0 = [Swp(id,7);Swp(id,6);Swp(id,5)];
    R = Quat2Rot(Swp(id,1:4)).';
    %R = [R(:,3),R(:,2),R(:,1)];
    
    H0 = [R,p0;0,0,0,1];
    H1 = [eye(3),p1;0,0,0,1];
    
    H = H0*H1;
    V(ii,:) = reshape(H(1:3,4),1,3);
end

mesh.Node = V;

end

%---------------------------------------------------------------- CURL MESH
function mesh = CurveMesh(mesh,Arg)

Node0 = mesh.Node;
sigma = max(Node0(:,3));
Ax = Arg{1};
dk = (Arg{2}*(pi/180)/sigma);
Celltmp = num2cell(Node0,2);

if strcmp(Ax,'x')
    CellDelta = cellfun(@(V) CurveCellOperation(V,dk,1e-12),...
    Celltmp,'UniformOutput',false);
elseif strcmp(Ax,'y')
    CellDelta = cellfun(@(V) CurveCellOperation(V,1e-12,dk),...
    Celltmp,'UniformOutput',false);
elseif strcmp(Ax,'z')
    
elseif strcmp(Ax,'PCC')
    dk1 = (Arg{2}*(pi/180)/sigma);
    dk2 = (Arg{3}*(pi/180)/sigma);
    CellDelta = cellfun(@(V) CurveCellOperation(V,dk1,dk2),...
    Celltmp,'UniformOutput',false);
elseif strcmp(Ax,'PCC+')
    dk1 = (Arg{2}*(pi/180))/(sigma*Arg{4});
    dk2 = (Arg{3}*(pi/180))/(sigma*Arg{4});
    Node0(:,3) = Arg{4}*Node0(:,3);
    Celltmp = num2cell(Node0,2);
    CellDelta = cellfun(@(V) CurveCellOperation(V,dk1,dk2),...
    Celltmp,'UniformOutput',false);
end

Node = vertcat(CellDelta{:});
mesh.Node = Node;

end

%---------------------------------------------------------------- CURL MESH
function mesh = TwistMesh(mesh,Arg)

Node0 = mesh.Node; 
Node = Node0;
Ax = Arg{1};
dk = (Arg{2}*(pi/180));

if strcmp(Ax,'x')
    Node(:,3) = Node0(:,1); Node(:,1) = Node0(:,3); 
    NodeC = [mean(Node(:,1)),mean(Node(:,2)),mean(Node(:,3))];
    Celltmp = num2cell(Node,2);
    upp = max(Node(:,3));
    low = min(Node(:,3));
    CellDelta = cellfun(@(V) TwistCellOperation(V,NodeC,dk,low,upp),...
    Celltmp,'UniformOutput',false);
    Node0 = vertcat(CellDelta{:});
    Node(:,3) = Node0(:,1); Node(:,1) = Node0(:,3); 
elseif strcmp(Ax,'y')
    Node(:,3) = Node0(:,2); Node(:,2) = Node0(:,3); 
    NodeC = [mean(Node(:,1)),mean(Node(:,2)),mean(Node(:,3))];
    Celltmp = num2cell(Node,2);
    upp = max(Node(:,3));
    low = min(Node(:,3));
    CellDelta = cellfun(@(V) TwistCellOperation(V,NodeC,dk,low,upp),...
    Celltmp,'UniformOutput',false);
    Node0 = vertcat(CellDelta{:});
    Node(:,3) = Node0(:,1); Node(:,1) = Node0(:,3);
elseif strcmp(Ax,'z')
    NodeC = [mean(Node(:,1)),mean(Node(:,2)),mean(Node(:,3))];
    Celltmp = num2cell(Node,2);
    upp = max(Node0(:,3));
    low = min(Node0(:,3));
    CellDelta = cellfun(@(V) TwistCellOperation(V,NodeC,dk,low,upp),...
    Celltmp,'UniformOutput',false);
    Node0 = vertcat(CellDelta{:});
    Node = Node0;
end

mesh.Node = Node;

end

%------------------------------------------------------- CURVATURE OPERATOR 
function V = CurveCellOperation(Node,kx,ky)

vx = Node(1);
vy = Node(2);
sigma = Node(3);

kappa = sqrt(kx^2 + ky^2);
theta = atan2(ky,kx); 

I3 = eye(3,3);
R = CurveRotationMatrix(sigma,kappa,theta);

p0 = (1/kappa)*[(1-cos(kappa*sigma))*cos(theta);
    (1-cos(kappa*sigma))*sin(theta);
    sin(kappa*sigma)];

p1 = [vx;vy;0];

H0 = [R,p0;0,0,0,1];
H1 = [I3,p1;0,0,0,1];

H = H0*H1;
V = reshape(H(1:3,4),1,3);

end

%------------------------------------------------ CURVATURE ROTATION MATRIX
function R = CurveRotationMatrix(sigma, kappa, theta)

alpha = kappa*sigma;

ca = cos(alpha);
sa = sin(alpha);
va = 1-ca;

kx = -sin(theta);
ky = cos(theta);

R = [(kx^2)*va + ca, kx*ky*va, ky*sa;
     kx*ky*va, (ky^2)*va + ca, -kx*sa;
     -ky*sa, kx*sa,  ca];

end

function V = TwistCellOperation(Node,NodeC,theta,sigmalow,sigmaupp)

vx = Node(1);
vy = Node(2);
vz = Node(3);
vxm = NodeC(1);
vym = NodeC(2);
sigma = (vz-sigmalow)/(sigmaupp - sigmalow);

I3 = eye(3,3);
R = TwistRotationMatrix(sigma, theta);

p0 = [vx-vxm;vy-vym;vz];
p1 = p0*0;
%p2 = [vxm;vym;0];

H0 = [R,p1;0,0,0,1];
H1 = [I3,p0;0,0,0,1];
%H2 = [I3,p2;0,0,0,1];

H = H0*H1;
V = reshape(H(1:3,4),1,3);

end

function R = TwistRotationMatrix(sigma, theta)
gamma = theta*sigma;

R = [ cos(gamma), sin(gamma), 0;
     -sin(gamma), cos(gamma), 0;
               0,          0, 1];
end


%---------------------------------------------------------------------- set
function R = Quat2Rot(q)
w = q(1); x = q(4); y = q(3); z = q(2);
Rxx = 1 - 2*(y^2 + z^2); Rxy = 2*(x*y - z*w); Rxz = 2*(x*z + y*w); 
Ryx = 2*(x*y + z*w); Ryy = 1 - 2*(x^2 + z^2); Ryz = 2*(y*z - x*w );
Rzx = 2*(x*z - y*w ); Rzy = 2*(y*z + x*w ); Rzz = 1 - 2 *(x^2 + y^2);

R = [Rxx, Rxy, Rxz; Ryx, Ryy, Ryz; Rzx, Rzy, Rzz];
end


