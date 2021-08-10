function Material = Elastosil28(D)
if nargin < 1, D = 10; end
Material=  YeohMaterial('C1',0.11,'C2',0.02,'C3',0,...
    'D1',D,'D2',1e-3,'D3',1e-3);


Material.Rho  = 1.132e-9;
Material.Zeta = 0.1;

end

