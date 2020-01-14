function d = dStrut(P,v1,v2,r)
a = r;
b = r;
c = sqrt(sum((v1 - v2).^2))/2;

dv = v2-v1; dv = dv/norm(dv);
R = PlanarProjection(dv);
P = transpose(R.'*P.');
% P = P + repmat([0,0,c],length(P),1);
% P = P + repmat(v1,length(P),1);
d = [abs(P(:,1))/a, abs(P(:,2))/b , abs(P(:,3))/c];
d = [d,max(d,[],2)-1+1e-6];
end


function R = PlanarProjection(n)
a = [0,0,1];
b = n(:)/norm(n);

v = cross(a,b); 
s = norm(v);
if s ~= 0
c = dot(a,b);
vs = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
R = eye(3,3) + vs + vs*vs*((1-c)/(s^2));
else
R = eye(3,3);
end

end
