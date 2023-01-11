function P = pMirror(P0,c) 
P = P0;
if norm(c(1)) ~= 0
    P(:,1) = abs(P(:,1));
elseif norm(c(2)) ~= 0
    P(:,2) = abs(P(:,2));
else
    P(:,3) = abs(P(:,3));
end
