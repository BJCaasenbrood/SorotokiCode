function [X,Y,Z,P] = MeshGridding(BdBox,Quality)
BdBox = BdBox*(1+1e-6);

if strcmp(Quality,'low'), N = 10;
elseif strcmp(Quality,'medium'), N = 50;
elseif strcmp(Quality,'high'), N = 100;
elseif isfloat(Quality), N = Quality;
else, N = 10;
end

x = linspace(BdBox(1),BdBox(2),N);
y = linspace(BdBox(3),BdBox(4),N);
if length(BdBox)>4, z = linspace(BdBox(5),BdBox(6),N);
else, z = 0; end

[X,Y,Z] = meshgrid(x,y,z); X = (X); Y = (Y); Z = (Z);
P = [X(:),Y(:),Z(:)]; 

end