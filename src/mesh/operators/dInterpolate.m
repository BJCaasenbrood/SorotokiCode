function d = dInterpolate(d1,d2,k)
% DINTERPOLATE the linear interpolation between two SDF fields. Usages:
%
%   a = 0.75;
%   d = dInterpolate(d1,d2,a)   %  takes 75% of d1, and 25% of d2.
%
%   Last edit: 15 July, 2022.

d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
d=[d,lerp(d1(:,end),d2(:,end),clamp(k,0,1))];