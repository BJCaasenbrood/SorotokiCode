function Material = Elastosil28(D)
if nargin < 1, D = 2; end
Material=  YeohMaterial('C1',0.11,'C2',0.02,'C3',0,...
    'D1',D,'D2',D,'D3',D);


Material.Rho  = 1.332e-9;
Material.Zeta = 0.2;
Material.Cfr  = 2e-6;

end

