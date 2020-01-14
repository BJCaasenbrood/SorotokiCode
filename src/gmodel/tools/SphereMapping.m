function UV = SphereMapping(r,a)
if nargin < 2, a = 1; end
m = 2*sqrt((r(:,1)).^2 + r(:,2).^2 + (abs(r(:,3)) + 1).^2);
m = repmat(m,[1,2]);
UV = a*(0.95./m).*r(:,1:2)+0.5;
UV = fliplr(UV);
end