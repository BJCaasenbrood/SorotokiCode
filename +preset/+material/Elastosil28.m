function Material = Elastosil28(D)
if nargin < 1, D = 1; end
Material=  Yeoh('C1',0.11,'C2',0.02,'C3',0,...
    'D1',D,'D2',D,'D3',D);

Material.params.Rho  = 1.332e-9;
Material.params.Zeta = 0.2;
end

