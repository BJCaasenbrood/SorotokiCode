function Material = Ecoflex0050(D)
if nargin < 1, D = 10; end
%Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    %'D1',D,'D2',10,'D3',10);
    
Material = Yeoh('C1',0.0145,'C2',0.0001,'C3',0,...
    'D1',D,'D2',20,'D3',30);

Material.params.Rho = 1070e-12;
Material.params.Zeta = 0.035;
Material.contact.TangentFriction = 0.75;
Material.contact.NormalReaction = 0.75;
end
