%------------------------------- GENERATE RANDOM POINTSET INSIDE HEMISPHERE 
function T = generateHemiSphereBase(Size,Bias,HalfBool)
if nargin < 3, HalfBool = true; end
T = zeros(Size,3);

for ii = 1:Size
    if HalfBool
        V1 = random(-1,1);
        V2 = random(-1,1);
        V3 = random(0,1);
        Kern = [V1,V2,V3];
    else
        V1 = random(-1,1);
        V2 = random(-1,1);
        V3 = random(-1,1);
        Kern = [V1,V2,V3];
    end
    Kern = Kern/norm(Kern);
    
    % power 1/2 is to ensure uniform sample in S3
    scale = lerp(Bias,1.0, (ii/Size)^0.5);
    T(ii,:) = Kern*scale;
end

end
