function plot3vector(p,n)
hold on; quiver3(p(:,1), p(:,2),p(:,3), n(:,1), n(:,2), n(:,3),'Color',col(2));
end
