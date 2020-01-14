function d = dCylinder(P,xc,yc,z1,z2,r)
d1 = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
d2 = z1-P(:,3);
d3 = P(:,3)-z2;

D1 = num2cell(d1); D2 = num2cell(d2); D3 = num2cell(d3);

D = cellfun(@(E,F,G) cases(E,F,G),D1,D2,D3,'UniformOutput',false);

d = [[d1,d2,d3], vertcat(D{:})];
end
%-------------------------------------------------------------------------%

function Dist = cases(d1,d2,d3)

if d1 > 0 && d2 >0
Dist = sqrt(d1^2 + d2^2);
elseif d1> 0 && d3>0
Dist = sqrt(d1^2 + d3^2);   
else
Dist = max([d1,d2,d3],[],2);  
end

end