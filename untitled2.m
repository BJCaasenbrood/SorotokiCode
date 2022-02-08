clr; beep off;
%% STL
obj = Gmodel('rubberducky.stl');
obj.Texture = diffuse(0.15);
FV = struct;
FV.faces = obj.Element;
FV.vertices = obj.NNode;

obj.bake.render();
view(30,30);

%% bounding box
BdBox = 1.25*round(boxhull(obj.Node)); BdBox(5) = BdBox(5) - 10;

x = linspace(BdBox(1),BdBox(2),150);
y = linspace(BdBox(3),BdBox(4),5);
z = linspace(BdBox(5),BdBox(6),150);
[X,Y,Z] = meshgrid(x,y,z);

V = obj.Node;
P = [X(:),Y(:),Z(:)];

%IDX = knnsearch(P,V);
d = DistancePointSet(P,V);
%d = sqrt(sum((P-V(IDX,:)).^2,2));
D = reshape(d,size(X));

%% slice
hold on;
h = slice(X,Y,Z,D,[],0,[]); h.LineStyle = 'none';
colormap(turbo());

function [dist,I] = DistancePointSet(PS1,PS2)

%d = zeros(size(PS2,1),size(PS1,1));

dist = zeros(size(PS2,1),1);
I    = zeros(size(PS2,1),1);

for el = 1:size(PS1,1)
    if size(PS1,2) == 3
        d = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2 + ...
            (PS1(el,3)-PS2(:,3)).^2 );
    else
        d = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    end
    
    [dist(el),I(el)] = min(d);
end



for el = 1:size(PS2,1)
    [dist(el),I(el)] = min(d(el,:));
    %[dist(el),I(el)] = min(d(el,:));
end
end