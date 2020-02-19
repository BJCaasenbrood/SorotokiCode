function Pc = Quads(BdBox,Nx,Ny)
if nargin < 3, Ny = Nx; end

B1 = BdBox(1); B2 = BdBox(2);
B3 = BdBox(3); B4 = BdBox(4);

DBx = (B2-B1)/Nx; DBy = (B4-B3)/Ny;
x = stepspace(B1+DBx/2,DBx,Nx);
y = stepspace(B3+DBy/2,DBy,Ny);

[X,Y] = meshgrid(x,y);

Pc = [X(:),Y(:)];
end