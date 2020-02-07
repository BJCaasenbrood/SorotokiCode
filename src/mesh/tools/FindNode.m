function x = FindNode(Node,varargin)

tol = BuildTolerance(Node); 
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
    case('Location');   x = Location(Node,varargin{2:end});
    case('SDF');        x = SignedFunction(Node,varargin{2:end},tol);
    case('Line');       x = Line(Node,varargin{2:end},tol);
    case('FloodFill');  x = FloodFill(Node,varargin{2:end});
end

end

function id = SignedFunction(Node,SDF,tol)
d = SDF(Node); d = d(:,end);
id = find( abs(d) < tol );
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
id = find(abs(Node(:,2) - BdBox(4))<tol);
end

function id = Bottom(Node,BdBox,tol)
id = find(abs(Node(:,2) - BdBox(3))<tol );
end

function id = Left(Node,BdBox,tol)
id = find(abs(Node(:,1) - BdBox(1))<tol);
end

function id = Right(Node,BdBox,tol)
id = find(abs(Node(:,1) - BdBox(2))<tol);
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

function eps = BuildTolerance(Node)
BdBox = BoundingBox(Node);
eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
end

function id = Location(Node,Node0,N)
if nargin < 3, N = 1; end
Loc = sqrt((Node(:,1)-Node0(1)).^2+(Node(:,2)-Node0(2)).^2);
[~,Loc] = sort(Loc); id = Loc(1:N);
end

function id = Line(Node,Line,tol)
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
list=zeros(n,1);
list(v0)=1;
neigh=(m(:,v0) & ~list);
stack = find(neigh);

while(~isempty(stack))
    v0 = stack(end);
    stack(end) = [];
    list(v0)=1;
    neigh=(m(:,v0) & ~list); 
    stack = unique( [stack; find(neigh)]);  
end
end

function BdBox = BoundingBox(Node)
BdBox = zeros(4,1);
BdBox(1) = min(Node(:,1));
BdBox(2) = max(Node(:,1));
BdBox(3) = min(Node(:,2));
BdBox(4) = max(Node(:,2));
end

