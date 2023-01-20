%SPOLYLINE - returns a signed distance function for a polyline
%
%   sdf = sPolyline(X) returns a signed distance function for a polyline defined
%   by a set of points X. Each row of X represents the x and y coordinates of a point 
%   on the polyline.
%
%   The output sdf is a SDF object that can be used to evaluate the signed
%   distance of a given point to the polyline.
%
%   Example:
%       X = linspace(pi, 3*pi, 100).';
%       Y = sin(X);
%       sdf = sPolyline([X(:), Y(:)]);
%       sdf.show();
%
%   See also SDF, SPOLYLINE, SLINE
%

function sdf = sPolyline(X)

% if size(X,2) == 2
%     xmax = 1.25*max(X(:,1));
%     xmin = 1.25*min(X(:,1));
%     ymax = 1.25*max(X(:,2));
%     ymin = 1.25*min(X(:,2));
% else
%     
% end

sdf = Sdf(@(P) sdfLine(P,X));
sdf.BdBox = [min(X(:,1)), max(X(:,1)),...
             min(X(:,2)), max(X(:,2))];

end
%------------------------------------------------------------- Vector Class
function d = sdfLine(P,X)

Tvec = TangentField(X);

if size(P,2) == 3
   P = P(:,[1,2]); 
end

[XY,D,I] = distance2curve(X,P,'linear');
s = linspace(0,1,size(X,1)).';
F = griddedInterpolant(s,Tvec);

N = F(I);

dr = XY - P;
dr = dr./sqrt((sum((dr.^2),2)));
dr(:,3) = 0; 
N(:,3)  = 0;

so  = cross(N,dr);
dir = zeros(size(so,1),1);

for ii = 1:size(so,1)
   dir(ii) = sign(-dot(so(ii,:), [0,0,1]));
   if dir(ii) == 0 && ii > 2
       dir(ii) = dir(ii - 1);
   end
end

d = dir.*[D,D];

end
%------------------------------------------------------------- Vector Class
function T = TangentField(PS1)
[~, Fy] = gradient(PS1);
T = [Fy(:,1), Fy(:,2)];
T = T ./ sqrt((sum((T.^2),2)));
end