%SCIRCLE - returns a signed distance function for a circle
%
%   sdf = sCircle(r) returns a signed distance function for a circle with 
%   radius r and centered at the origin (0, 0).
%
%   sdf = sCircle(r, [x, y]) returns a signed distance function for a circle 
%   with radius r and centered at point (x, y).
%
%   The output sdf is a SDF object that can be used to evaluate the signed
%   distance of a given point to the circle.
%
%   Example:
%       sdf = sCircle(1, [0, 0])
%       point = [-0.5, 0.5];
%       dist = sdf.eval(point)
%
%   See also SDF, DCIRCLE

function sdf = sCircle(varargin)
r  = 1;
xc = 0;
yc = 0;

if nargin == 1
    r  = abs(varargin{1}); 
elseif nargin > 1
    r = abs(varargin{1});
    if numel(varargin) == 2
        v  = varargin{2};
        xc = v(1); 
        yc = v(2);
    else
        error('Translation should be of size 2 in sCircle(r,[x,y])');
    end    
end

eps = 1e-4*r;
sdf = Sdf(@(P) dCircle(P,xc,yc,r));
sdf.BdBox = [xc - r - eps, xc + r + eps, ...
             yc - r - eps, yc + r + eps];

end

function d = dCircle(P,xc,yc,r)
d = sqrt((P(:,1) - xc).^2 + (P(:,2) - yc).^2) - r;
end