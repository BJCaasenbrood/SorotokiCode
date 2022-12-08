function plot3vector(p,n)
% PLOT3VECTOR plots an 3D-vector. Usage:
%    pos = [1,2,3];
%    dir = [1,0,1];
%    plotvector3(pos,dir);
%
%   Last edit: Dec 3, 2022
hold on; 
quiver3(p(:,1), p(:,2),p(:,3), ...
        n(:,1), n(:,2), n(:,3),'Color',col(2));
end
