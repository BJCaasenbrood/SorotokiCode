function d = dDiff(d1,d2) % max(d1,-d2)
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
d=[d,max(d1(:,end),-d2(:,end))];