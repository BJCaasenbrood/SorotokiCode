function Material = Dragonskin30(D)
if nargin < 1, D = 1; end

Material = YeohMaterial('C1',0.1536,'C2',0.0016,'C3',0,...
    'D1',D,'D2',10,'D3',10);

%Material = YeohMaterial('C1',1.19e-3,'C2',2.3028e-2,'C3',1e-12,...
%    'D1',D,'D2',1e-3,'D3',1e-3);
end

