%SRECTANGLE creates a signed distance function (SDF) of a rectangle
%
%   sdf = sRectangle() creates a default rectangle with its bottom-left corner 
%   at (0,0) and its top-right corner at (1,1)
% 
%   sdf = sRectangle(L) creates a square with side length L
% 
%   sdf = sRectangle(W, H) creates a rectangle with width W and height H
% 
%   sdf = sRectangle([x1,y1], [x2,y2]) creates a rectangle with bottom-left 
%   corner at (x1, y1) and top-right corner at (x2, y2)
% 
%   sdf = sRectangle(x1, x2, y1, y2) creates a rectangle with bottom-left
%   corner at (x1, y1) and top-right corner at (x2, y2)
%
% The function outputs a 'sdf' object, which is a representation of the SDF of the rectangle. 
%
%   See also SCIRCLE, SCUBE

function sdf = sRectangle(varargin)
x1 = 0;
x2 = 1;
y1 = 0;
y2 = 1;

if nargin == 1
   x2 = varargin{1};
   y2 = varargin{1};
elseif nargin == 2
   if numel(varargin{1}) == 1
       x2 = varargin{1};
       y2 = varargin{2};
   elseif numel(varargin{1}) == 2
       p0 = varargin{1}; p1 = varargin{2};
       x1 = p0(1);
       y1 = p0(2);
       x2 = p1(1);
       y2 = p1(2);
       
       if x1 == x2
           error('Rectangle has zero width');
       elseif y1 == y2
           error('Rectangle has zero height');
       end
   else
      error('Wrong input for sRectangle');
   end
elseif nargin == 4
   x1 = varargin{1};
   y1 = varargin{3};
   x2 = varargin{2};
   y2 = varargin{4};
end

eps = 1e-4 * norm([x1; x2; y1; y2]);
sdf = Sdf(@(P) sdfRectangle(P,x1,x2,y1,y2));
sdf.BdBox = [x1 - eps, x2 + eps, y1 - eps, y2 + eps];
end

function d = sdfRectangle(P,x1,x2,y1,y2)
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2];
d = [d, max(d,[],2)];
end