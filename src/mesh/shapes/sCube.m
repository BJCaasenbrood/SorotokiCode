function sdf = sCube(x1,x2,y1,y2,z1,z2)
if nargin < 2
    if numel(x1) == 1
        a = x1;
        x1 = -a;
        x2 = a; 
        y1 = -a;
        y2 = a;
        z1 = -a;
        z2 = a;
    else
        X = x1;
        x1 = X(1);
        x2 = X(2);
        y1 = X(3);
        y2 = X(4);
        z1 = X(5);
        z2 = X(6);
    end
end

sdf = Sdf(@(P) sdfCube(P,x1,x2,y1,y2,z1,z2));
sdf.BdBox = [x1-1e-6,x2+1e-6,....
             y1-1e-6,y2+1e-6,...
             z1-1e-6,z2+1e-6];
end

function d = sdfCube(P,x1,x2,y1,y2,z1,z2)
d = [x1-P(:,1) + 1e-12, P(:,1)-x2, y1-P(:,2)+1e-12, P(:,2)-y2, z1-P(:,3)+1e-12, P(:,3)-z2];
d = [d,max(d,[],2)];
end
