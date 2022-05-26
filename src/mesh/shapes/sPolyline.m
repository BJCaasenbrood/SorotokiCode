function sdf = sPolyline(X)

if size(X,2) == 2
    xmax = max(X(:,1));
    xmin = min(X(:,1));
    ymax = max(X(:,2));
    ymin = min(X(:,2));
else
    
end

sdf = Sdf(@(P) sdfLine(P,X));
sdf.BdBox = [xmin,xmax,ymin,ymax];

n = size(X,1);
T = TangentField(X);

sdf.Node    = X + 1e-3*T;
sdf.Element = [(1:n-1).',(2:n).'];

end
%------------------------------------------------------------- Vector Class
function d = sdfLine(P,X)

Tvec = TangentField(X);
%[D,I] = DistancePointSet(X,P);
[XY,D,I] = distance2curve(X,P,'linear');
s = linspace(0,1,size(X,1)).';
F = griddedInterpolant(s,Tvec);

N = F(I);

%N  = Tvec(I,:);
dr = XY - P;
dr = dr./sqrt((sum((dr.^2),2)));
dr(:,3) = 0; N(:,3)  = 0;

so  = cross(N,dr);
dir = zeros(length(so),1);
for ii = 1:length(so)
   dir(ii) = sign(-dot(so(ii,:),[0,0,1]));
   if dir(ii) == 0
       dir(ii) = dir(ii-1);
   end
end

d = dir.*[D,D];

end
%------------------------------------------------------------- Vector Class
function [dist,I] = DistancePointSet(PS1,PS2)
d = zeros(size(PS2,1),size(PS1,1));
for el = 1:size(PS1,1)
    if size(PS1,2) == 3
        dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2 + ...
            (PS1(el,3)-PS2(:,3)).^2 );
    else
        dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    end
    
    d(:,el) = dist;
end

dist = zeros(size(PS2,1),1);
I = zeros(size(PS2,1),1);
for el = 1:size(PS2,1)
    [dist(el),I(el)] = min(d(el,:));
    %[dist(el),I(el)] = min(d(el,:));
end
end
%------------------------------------------------------------- Vector Class
function T = TangentField(PS1)
%FF = [diff(PS1(:,1)),diff(PS1(:,2))];
%Fy = [FF;FF(end,:)];
[~, Fy] = gradient(PS1);
T = [Fy(:,1),Fy(:,2)];
T = T./sqrt((sum((T.^2),2)));
end