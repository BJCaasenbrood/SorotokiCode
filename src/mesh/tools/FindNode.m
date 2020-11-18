function x = FindNode(Node,varargin)

tol   = BuildTolerance(Node); 
BdBox = BoundingBox(Node);

Request = varargin{1}; 

switch(Request)
    case('SW');         x = SouthWest(Node,BdBox,tol);
    case('SE');         x = SouthEast(Node,BdBox,tol);
    case('NW');         x = NorthWest(Node,BdBox,tol);
    case('NE');         x = NorthEast(Node,BdBox,tol);
    case('Top');        x = Top(Node,BdBox,tol);
    case('Bottom');     x = Bottom(Node,BdBox,tol);
    case('Left');       x = Left(Node,BdBox,tol);
    case('Right');      x = Right(Node,BdBox,tol);
    case('Middle');     x = Middle(Node,BdBox,varargin{end});
    case('TopMid');     x = TopMid(Node,BdBox);
    case('BottomMid');  x = BottomMid(Node,BdBox);
    case('LeftMid');    x = LeftMid(Node,BdBox,varargin{end});
    case('RightMid');   x = RightMid(Node,BdBox);
    case('TopHalf');    x = TopHalf(Node,BdBox,tol);
    case('BotHalf');    x = BotHalf(Node,BdBox,tol);
    case('Location');   x = Location(Node,BdBox,varargin{2:end});
    case('SDFE');       x = SignedFunctionEdge(Node,varargin{2:end},tol);
    case('SDF');        x = SignedFunction(Node,varargin{2:end},tol);
    case('Line');       x = Line(Node,varargin{2:end},tol);
    case('Box');        x = BoxSelect(Node,varargin{2:end},tol);
    case('FloodFill');  x = FloodFill(Node,varargin{2:end});
end

end

function id = SignedFunctionEdge(Node,SDF,tol)
d = SDF(Node); d = d(:,end);
id = find( abs(d) < tol );
end


function id = SignedFunction(Node,SDF,tol)
d = SDF(Node); d = d(:,end);
id = find( (d) < 0.01*tol );
end

function id = SouthEast(Node,BdBox,tol)
id = find(abs(Node(:,1)-BdBox(2)) < tol & abs(Node(:,2)-BdBox(3)) < tol );
end

function id = SouthWest(Node,BdBox,tol)
id = find(abs(Node(:,1)-BdBox(1)) < tol & abs(Node(:,2)) < tol );
end

function id = NorthWest(Node,BdBox,tol)
id = find(abs(Node(:,1)-BdBox(1)) < tol & abs(Node(:,2)-BdBox(4)) < tol );
end

function id = NorthEast(Node,BdBox,tol)
id = find(abs(Node(:,1)-BdBox(2)) < tol & abs(Node(:,2)-BdBox(4)) < tol );
end

function id = Top(Node,BdBox,tol)
if length(BdBox) == 4, id = find(abs(Node(:,2) - BdBox(4))<tol);
else, id = find(abs(Node(:,3) - BdBox(6))<tol);
end
end

function id = Bottom(Node,BdBox,tol)
if length(BdBox) == 4, id = find(abs(Node(:,2) - BdBox(3))<tol );
else, id = find(abs(Node(:,3) - BdBox(5))<tol);
end
end

function id = Left(Node,BdBox,tol)
if length(BdBox) == 4, id = find(abs(Node(:,1) - BdBox(1))<tol);
else, id = find(abs(Node(:,1) - BdBox(1))<tol);
end
end

function id = Right(Node,BdBox,tol)
if length(BdBox) == 4, id = find(abs(Node(:,1) - BdBox(2))<tol);
else, id = find(abs(Node(:,1) - BdBox(2))<tol);
end
end

function id = TopHalf(Node,BdBox,tol)
id = find(Node(:,2) - ((BdBox(4)+BdBox(3))/2 + BdBox(3))>tol);
end

function id = BotHalf(Node,BdBox,tol)
id = find(Node(:,2) - ((BdBox(4)+BdBox(3))/2 + BdBox(3))<tol);
end

function id = Middle(Node,BdBox,N)
Mid = sqrt((Node(:,1)-(BdBox(1)+BdBox(2))*0.5).^2+...
    (Node(:,2)-(BdBox(3)+BdBox(4))*0.5).^2);
[~,Mid] = sort(Mid); id = Mid(1:N);
end

function id = RightMid(Node,BdBox)
Mid = sqrt((Node(:,1)-BdBox(2)).^2+...
    (Node(:,2)-(BdBox(3)+BdBox(4))*0.5).^2);
[~,Mid] = sort(Mid); id = Mid(1);
end

function id = LeftMid(Node,BdBox,N)
if nargin < 3 || ischar(N), N = 1; end
Mid = sqrt((Node(:,1)-BdBox(1)).^2+...
    (Node(:,2)-(BdBox(3)+BdBox(4))*0.5).^2);
[~,Mid] = sort(Mid); id = Mid(1:N);
end

function id = TopMid(Node,BdBox,N)
if nargin < 3, N = 1; end
Mid = sqrt((Node(:,1)-(BdBox(1)+BdBox(2))*0.5).^2+...
    (Node(:,2)-BdBox(4)).^2);
[~,Mid] = sort(Mid); id = Mid(1:N);
end

function id = BottomMid(Node,BdBox,N)
if nargin < 3, N = 1; end
Mid = sqrt((Node(:,1)-(BdBox(1)+BdBox(2))*0.5).^2+...
    (Node(:,2)-BdBox(3)).^2);
[~,Mid] = sort(Mid); id = Mid(1:N);
end


function eps = BuildTolerance(Node)
BdBox = BoundingBox(Node);
eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
end

function id = Location(Node,BdBox,Node0,N)
if nargin < 4, N = 1; end
if length(BdBox) == 4
    Loc = sqrt((Node(:,1)-Node0(1)).^2+(Node(:,2)-Node0(2)).^2);
elseif length(BdBox) == 6
    Loc = sqrt((Node(:,1)-Node0(1)).^2+ ...
        (Node(:,2)-Node0(2)).^2 + (Node(:,3)-Node0(3)).^2);
end

[~,Loc] = sort(Loc); id = Loc(1:N);
end

function id = Line(Node,Line,tol)
d = dLine(Node,Line(1)-eps,Line(2)+eps,Line(3)-eps,Line(4)+eps);
id = find(abs(d(:,end))<tol);
end

function id = BoxSelect(Node,Line,tol)
d = dRectangle(Node,Line(1)-eps,Line(2)+eps,Line(3)-eps,Line(4)+eps);
id = find(abs(d(:,end))<tol);
end


function id = FloodFill(Node,Fem,rho)

m = Fem.Mesh.get('Adjecency');

blclist = find(rho>Fem.get('VoidTolerance'));
m(blclist,:) = 0;
m(:,blclist) = 0;

v0 = Fem.get('PressureCell');

x = FloodFillOne(m,v0(1));
id = 1:length(Node);
id = id(logical(x));

end

function list = FloodFillOne(m,v0)

[n,~]=size(m);
L = 1:n;
list=zeros(n,1);
list(v0)=1;
neigh=m(:,v0).*(~list);%(m(:,v0) & ~list);
stack = find(neigh);

while(~isempty(stack))
    v0 = stack(end);
    stack(end) = [];
    list(v0)=1;
    neigh= m(:,v0).*(~list);
     %neigh= (m(:,v0) & ~list);
    stack = sort([stack; find(neigh)]); 
    stack = stack(diff(stack(:))>0);
end

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

