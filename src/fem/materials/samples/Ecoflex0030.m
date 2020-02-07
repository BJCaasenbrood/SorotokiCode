function Material = Ecoflex0030(Model)
if nargin < 1
Material=  YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',10,'D2',10,'D3',10);
else
if strcmp(Model,'Yeoh')
Material=  YeohMaterial('C1',7.61e-3,'C2',2.42e-4,'C3',-6.2e-7,...
     'D1',1.0,'D2',1.0,'D3',1.0);
%Material = YeohMaterial('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
%     'D1',15,'D2',10,'D3',10);
elseif strcmp(Model,'Mooney')
end
end

end

