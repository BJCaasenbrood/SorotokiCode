function [Normal, Center] = TriangleNormal(Node,Element)
N = zeros(length(Element),3);
C = zeros(length(Element),3);

for i = 1:length(Element)
    v1 = Node(Element(i,2),:) - Node(Element(i,1),:);
    v2 = Node(Element(i,3),:) - Node(Element(i,1),:);
    
    tmp = cross(v1,v2); 
    tmp = tmp/sqrt(tmp(1)^2 + tmp(2)^2 + tmp(3)^2);
    for j = 1:3
       if isnan(tmp(1)), tmp(1) = 0; end
       if isnan(tmp(2)), tmp(2) = 0; end
       if isnan(tmp(3)), tmp(3) = 1; end
    end
    
    N(i,:) = tmp;

    p1 = Node(Element(i,1),:);
    p2 = Node(Element(i,2),:);
    p3 = Node(Element(i,3),:);
    
    C(i,:) = (1/3).*(p1 + p2 + p3);
end

Normal = N;
Center = C;
end
