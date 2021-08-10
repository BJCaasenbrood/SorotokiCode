%------------------------------------------ PLANAR PROJECTION OF 3D-POLYGON
function R = PlanarProjection(n)
b = n(:)/norm(n); a = [0;0;1];

v = cross(a,b); 
s = norm(v);
c = dot(a,b);
vs = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
R = eye(3,3) + vs + vs*vs*((1-c)/(s^2));
end