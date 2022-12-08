function Node = depth2nodes(pointcloud,depth,varargin)

if isempty(varargin)
    BdBox = [-1e6,1e6,-1e6,1e6,-1e6,1e6];
else
    BdBox = varargin{1};
end

%sdf = sCube(BdBox);

if (depth.logical())% && color.logical())

    %pointcloud.map_to(color);
    points = pointcloud.calculate(depth);

    % Adjust frame CS to matlab CS
    vertices = points.get_vertices();
    Y = -vertices(:,1,1);
    Z = -vertices(:,2,1);
    X = vertices(:,3,1);
end

P = [X(:),Y(:),Z(:)];
x1 = BdBox(1); x2 = BdBox(2);
y1 = BdBox(3); y2 = BdBox(4);
z1 = BdBox(5); z2 = BdBox(6);
d = [x1-P(:,1) + 1e-12, P(:,1)-x2, y1-P(:,2)+1e-12, P(:,2)-y2, z1-P(:,3)+1e-12, P(:,3)-z2];
d = max(d,[],2);
Y = d<0;
%Y = sdf.intersect(Node0);
%sum(Y);
Node = P(Y,:);

end

