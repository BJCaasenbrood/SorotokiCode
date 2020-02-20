function Material = Elastosil
d = 1.0;
Material=  YeohMaterial('C1',0.11,'C2',0.02,'C3',0,...
    'D1',d,'D2',1e3,'D3',1e3);
end

