%------------------------------------------------- CHECK UV-TEXTURE MAPPING
function [x,y] = FitMaterialCap(x,y,Nx,Ny)
%#codegen
R = floor(mean([Nx,Ny])/2);
Xm = floor(Nx/2); Ym = floor(Ny/2);
Delta = [x - Xm, y - Ym]; 
Theta = atan2(Delta(1),Delta(2));
dR = round(sqrt(Delta(1)^2 +  Delta(2)^2));

if dR >= R
   x = floor(Xm + sin(Theta)*0.99*R);
   y = floor(Ym + cos(Theta)*0.99*R);
end

if real(x) >= Nx, x = Nx-2; end
if real(y) >= Ny, y = Ny-2; end
if real(x) <= 1, x = 2; end
if real(y) <= 1, y = 2; end

if isnan(x), x = round(0.5*Nx)-5; end
if isnan(y), y = round(0.5*Ny)-5; end

x = real(x); y = real(y);
end