function Material = Dragonskin10(D)
if nargin < 1, D = 1; end
% Material = YeohMaterial('C1',36e-3,'C2',0.25e-3,'C3',0.023e-3,...
%     'D1',D,'D2',1e-3,'D3',1e-3);
Material = Yeoh('C1',36e-3,'C2',0.25e-3,'C3',0.023e-3,...
     'D1',D,'D2',D,'D3',D);
% Material = YeohMaterial('C1',0.0355,'C2',0.0110,'C3',0.0040,...
%      'D1',D,'D2',D,'D3',D);
% Material = YeohMaterial('C1',0.0245,'C2',-0.0001,'C3',0.0,...
%       'D1',D,'D2',D,'D3',D);

Material.params.Rho = 1070e-12;
Material.params.Zeta = 0.1;
end

