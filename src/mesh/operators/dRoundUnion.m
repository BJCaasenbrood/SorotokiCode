function d = dRoundUnion(d1,d2,k) % min(d1,d2)
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];

f1 = d1(:,end); F1 = num2cell(f1); 
f2 = d2(:,end); F2 = num2cell(f2);

dblend = cellfun(@(E,V) Roundoperation(E,V,k),F1,F2,'UniformOutput',false);

d=[d,vertcat(dblend{:})];

end

function d = Roundoperation(f1,f2,k)
intspace = min([f1-k,f2-k],0);

% h = max([k-abs(f1-f2), 0.0 ]);
% d = min([ f1, f2 ]) - h*h*h/(6.0*k*k);
end