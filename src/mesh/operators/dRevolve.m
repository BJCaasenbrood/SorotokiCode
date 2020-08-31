function y = dRevolve(x,rot,flip)

if nargin == 3
    if islogical(flip), x = (roty(pi/2)*x.').';
    else, x = (roty(flip)*x.').';
    end
end

%y = x - repmat(move(:)',length(x),1);
[th,r,z] = cart2pol(x(:,1),x(:,2),x(:,3));
if nargin == 1, y = [r(:),z(:)];
else
[X,Y,Z] = pol2cart(clamp(abs(th),0,rot),r,z);
[~,r,z] = cart2pol(X,Y,Z);    
y = [r(:),z(:)];
end

end

