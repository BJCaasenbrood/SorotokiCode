%------------------------------- GENERATE RANDOM POINTSET INSIDE HEMISPHERE 
function T = GenerateHemiSphereBase(Size,Bias,HalfBool)
if nargin < 3, HalfBool = true; end
T = zeros(Size,3);

for ii = 1:Size
    if HalfBool
    Kern = [random(-1,1), random(-1,1), random(0,1)];
    else
    Kern = [random(-1,1), random(-1,1), random(-1,1)];
    end
    Kern = Kern/norm(Kern);
    scale = lerp(Bias,1.0, (ii/Size)^2);
    T(ii,:) = Kern*scale;
end

end
