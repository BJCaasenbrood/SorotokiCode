%SCAP - returns a 2D signed distance function for a cupped rectangle
%
%   sdf = sCup(r, x) creates a 2D signed distance function for a cupped circle 
%   with radius r and height x. If x is negative, the shape is reversed.
%
%   Example: 
%       sdf = sCup(1, 2)
%    
%   See also SDF, SCAP, SCIRCLE

function sdf = sCup(varargin)
r  = 1;
xc = 0;
yc = 0;
x  = 1;

if nargin == 1
    r = varargin{1}; 
    x = r;
elseif nargin == 2
    r = varargin{1}; 
    x = varargin{2}; 
end

eps = 1e-4*r;

if x > 0
    sdf = Sdf(@(P) dDiff(dRectangle(P,xc-r,xc+r,yc-x,yc),...
        dCircle(P,xc,yc,r)));

    sdf.BdBox = [xc-r-eps,xc+r+eps,yc-x-eps,yc+eps];
else
    sdf = Sdf(@(P) dDiff(dRectangle(P,xc-r,xc+r,yc,yc+abs(x)),...
        dCircle(P,xc,yc,r)));

    sdf.BdBox = [xc-r-eps,xc+r+eps,yc-eps,yc+abs(x)+eps];
end

end

function d = dCircle(P,xc,yc,r)
    d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
    %d = [d,d];
end

function d = dRectangle(P,x1,x2,y1,y2)
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2];
d = max(d,[],2);
end