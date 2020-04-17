function Material = Ecoflex0030(D)
if nargin < 1, D = 1.0; end
Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',D,'D2',20,'D3',30);
end

