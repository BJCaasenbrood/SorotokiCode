function Material = Ecoflex0050(D)
if nargin < 1, D = 3; end
%Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    %'D1',D,'D2',10,'D3',10);
    
Material = YeohMaterial('C1',0.0145,'C2',0.0001,'C3',0,...
    'D1',D,'D2',10,'D3',10);
end

