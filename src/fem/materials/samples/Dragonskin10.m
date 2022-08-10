function Material = Dragonskin10(D)
if nargin < 1, D = 3; end
% Material = YeohMaterial('C1',36e-3,'C2',0.25e-3,'C3',0.023e-3,...
%     'D1',D,'D2',1e-3,'D3',1e-3);
Material = YeohMaterial('C1',36e-3,'C2',0.25e-3,'C3',0.023e-3,...
    'D1',D,'D2',2*D,'D3',3*D);

Material.Rho  = 1.07e-9;
Material.Zeta = 1;
Material.Cfr  = 1e-5;

end

