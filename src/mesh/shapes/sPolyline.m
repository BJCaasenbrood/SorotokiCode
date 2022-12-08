function sdf = sPolyline(X)

if size(X,2) == 2
    xmax = 1.25*max(X(:,1));
    xmin = 1.25*min(X(:,1));
    ymax = 1.25*max(X(:,2));
    ymin = 1.25*min(X(:,2));
else
    
end

sdf = Sdf(@(P) sdfLine(P,X));
sdf.BdBox = [xmin,xmax,ymin,ymax];

n = size(X,1);
T = TangentField(X);

sdf.Node    = X + 1e-3*T;
sdf.Element = [(1:n-1).',(2:n).'];
sdf.Center  = sdf.centerofmass;

end
%------------------------------------------------------------- Vector Class
function d = sdfLine(P,X)

Tvec = TangentField(X);
%[D,I] = DistancePointSet(X,P);

if size(P,2) == 3,
   P = P(:,[1,2]); 
end

[XY,D,I] = distance2curve(X,P,'linear');
s = linspace(0,1,size(X,1)).';
F = griddedInterpolant(s,Tvec);

N = F(I);

%N  = Tvec(I,:);
dr = XY - P;
dr = dr./sqrt((sum((dr.^2),2)));
dr(:,3) = 0; 
N(:,3)  = 0;

so  = cross(N,dr);
dir = zeros(size(so,1),1);
for ii = 1:size(so,1)
   dir(ii) = sign(-dot(so(ii,:),[0,0,1]));
   if dir(ii) == 0 && ii > 2
       dir(ii) = dir(ii-1);
%    else
%        dir(ii) = 0;
   end
end

d = dir.*[D,D];

end
%------------------------------------------------------------- Vector Class
function T = TangentField(PS1)
%FF = [diff(PS1(:,1)),diff(PS1(:,2))];
%Fy = [FF;FF(end,:)];
[~, Fy] = gradient(PS1);
T = [Fy(:,1),Fy(:,2)];
T = T./sqrt((sum((T.^2),2)));
end