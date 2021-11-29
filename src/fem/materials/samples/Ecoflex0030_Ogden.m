function Material = Ecoflex0030_Ogden

Material = OgdenMaterial('C1',1.241,'C2',7.89e-9,'A1',3.034,...
    'A2',13.02,'Bulk',1e36);

Material.Rho  = 1.07e-9;
Material.Zeta = 0.1;

%Material = YeohMaterial('C1',0.0084,'C2',0.0001,'C3',0,...
%    'D1',D,'D2',20,'D3',30);
end

