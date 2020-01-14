function plotvector(p,n)
n = n*norm(p)*0.3; % end position of normal vector

%quiver3 syntax: quiver3(x,y,z,u,v,w)
quiver3(p(1), p(2), p(3), n(1), n(2), n(3),'Color',ColorScheme(2));
hold on; plot3(p(1),p(2),p(3),'y.','markersize',10);
end
