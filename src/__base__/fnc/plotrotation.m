function plotrotation(p,R,varargin)
%n = n*norm(p)*0.3; % end position of normal vector

s = mean(abs(axis));
R = 1.25*R*s;

if numel(p) == 3
%quiver3 syntax: quiver3(x,y,z,u,v,w)
hold on;
for ii = 1:3
quiver3(p(1), p(2), p(3), R(ii,1), R(ii,2), R(ii,3),varargin{:});
end
%hold on; plot3(p(1),p(2),p(3),'y.','markersize',10);
elseif numel(p) == 2
hold on; 
for ii = 1:2
quiver(p(1), p(2), R(1,ii), R(2,ii),'k',varargin{:});
end

elseif size(p,2) == 3
hold on; handle = quiver3(p(:,1), p(:,2), p(:,3), ...
    n(:,1), n(:,2), n(:,3), varargin{:});
end
end
