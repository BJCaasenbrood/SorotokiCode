function d = dUnionBlend(d1,d2,k) % min(d1,d2)
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];

f1 = d1(:,end); F1 = num2cell(f1); 
f2 = d2(:,end); F2 = num2cell(f2);

dblend = cellfun(@(E,V) BlendOperation(E,V,k),F1,F2,'UniformOutput',false);

d=[d,vertcat(dblend{:})];
end

%-------------------------------------------------------------------------%
function d = BlendOperation(f1,f2,k)
h = max( k-abs(f1-f2), 0 );
d = max( f1, f2 ) - (h^3)/(6.0*k^2);
end