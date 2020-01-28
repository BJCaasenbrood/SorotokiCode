function plotvector(p,n)
%n = n*norm(p)*0.3; % end position of normal vector

if length(p) == 3
%quiver3 syntax: quiver3(x,y,z,u,v,w)
quiver3(p(1), p(2), p(3), n(1), n(2), n(3),'Color','r');
%hold on; plot3(p(1),p(2),p(3),'y.','markersize',10);
else
hold on; quiver(p(1), p(2), n(1), n(2),'Color','r');
end
end
