function d = dInterpolate(d1,d2,k)
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
d=[d,lerp(d1(:,end),d2(:,end),clamp(k,0,1))];