function d = dDiff(d1,d2) 
%DDIFF the difference operator between two SDF fields which transelates to
%   dDiff = max(d1,-d2).
%
%   Last edit: July 15, 2022.

d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
d=[d,max(d1(:,end),-d2(:,end))];