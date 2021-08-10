function d = dSmoothUnion(d1,d2,k) % min(d1,d2)
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];

f1 = d1(:,end); F1 = num2cell(f1); 
f2 = d2(:,end); F2 = num2cell(f2);

dblend = cellfun(@(E,V) Smoothoperation(E,V,k),F1,F2,'UniformOutput',false);

d=[d,vertcat(dblend{:})];

end

function d = Smoothoperation(f1,f2,k)
h = max([k-abs(f1-f2), 0.0 ]);
d = min([ f1, f2 ]) - h*h*h/(6.0*k*k);
end