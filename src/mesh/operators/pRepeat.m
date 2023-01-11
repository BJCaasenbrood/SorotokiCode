function P = pRepeat(P0,c) 
if numel(c) == 3
P = zeros(size(P0,1),3);
    P(:,1) = mod(P0(:,1),c(1));
    P(:,2) = mod(P0(:,2),c(2));
    P(:,3) = mod(P0(:,3),c(3));
else
    P = zeros(size(P0,1),2);
    P(:,1) = mod(P0(:,1),c(1));
    P(:,2) = mod(P0(:,2),c(2));
end
end