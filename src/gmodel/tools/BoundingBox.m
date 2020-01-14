function A = BoundingBox(Node)
Vx = Node(:,1); Vy = Node(:,2); Vz = Node(:,3);

A = [floor(min(Vx)), ceil(max(Vx)),...
     floor(min(Vy)), ceil(max(Vy)),...
     floor(min(Vz)), ceil(max(Vz))];

end

