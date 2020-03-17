function plotpoint(p,varargin)
%n = n*norm(p)*0.3; % end position of normal vector

if length(p) == 3
%quiver3 syntax: quiver3(x,y,z,u,v,w)
plot3(p(1), p(2), p(3),varargin{:});
%hold on; plot3(p(1),p(2),p(3),'y.','markersize',10);
else
hold on; plot3(p(:,1), p(:,2),varargin{:});
end
end
