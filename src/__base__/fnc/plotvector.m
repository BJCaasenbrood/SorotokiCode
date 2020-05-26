function handle = plotvector(p,n,varargin)
%n = n*norm(p)*0.3; % end position of normal vector

if numel(p) == 3
%quiver3 syntax: quiver3(x,y,z,u,v,w)
hold on;
handle = quiver3(p(1), p(2), p(3), n(1), n(2), n(3),varargin{:});
%hold on; plot3(p(1),p(2),p(3),'y.','markersize',10);
elseif size(p,2) == 2
hold on; handle = quiver(p(:,1), p(:,2), n(:,1), n(:,2),varargin{:});
elseif size(p,2) == 3
hold on; handle = quiver3(p(:,1), p(:,2), p(:,3), ...
    n(:,1), n(:,2), n(:,3), varargin{:});
end
end
