function Material = Ecoflex0030(D)
if nargin < 1, D = 2; end
Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',D,'D2',D,'D3',D);

Material.Rho  = 1070e-12;
Material.Zeta = 0.03;
Material.Cfr  = 5e-6;
end

