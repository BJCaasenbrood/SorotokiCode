function Pc = Hexahedron(BdBox,Nx,Ny,Nz)
if nargin < 3, Ny = Nx; Nz = Nx; 
csqrt = @(x) x^(1/3);
Nx = ceil(csqrt(Nx));
Ny = ceil(csqrt(Ny));
Nz = ceil(csqrt(Nz));
end

B1 = BdBox(1); B2 = BdBox(2);
B3 = BdBox(3); B4 = BdBox(4);
B5 = BdBox(5); B6 = BdBox(6);

DBx = (B2-B1)/Nx; DBy = (B4-B3)/Ny; DBz = (B6-B5)/Nz;
x = stepspace(B1+DBx/2,DBx,Nx);
y = stepspace(B3+DBy/2,DBy,Ny);
z = stepspace(B5+DBz/2,DBz,Nz);

[X,Y,Z] = meshgrid(x,y,z);

Pc = [X(:),Y(:),Z(:)];
end