function sdf = sCube(x1,x2,y1,y2,z1,z2)
if nargin < 2
    if numel(x1) == 1
        a = x1; x1 = -a; x2 = a; 
        y1 = -a; y2 = a; z1 = -a;
        z2 = a;
    elseif numel(x1) == 3
        X = x1; x1 = 0; x2 = X(1);
        y1 = 0; y2 = X(2);
        z1 = 0; z2 = X(3);
    else
        X = x1; x1 = X(1); x2 = X(2);
        y1 = X(3); y2 = X(4);
        z1 = X(5); z2 = X(6);
    end
end

sdf = Sdf(@(P) sdfCube(P,x1,x2,y1,y2,z1,z2));
sdf.BdBox = [x1-1e-6,x2+1e-6,....
             y1-1e-6,y2+1e-6,...
             z1-1e-6,z2+1e-6];
         
sdf.Node = [x1 x2 x2 x1 x1 x2 x2 x1;
            y1 y1 y2 y2 y1 y1 y2 y2;
            z1 z1 z1 z1 z2 z2 z2 z2].';
   
sdf.Element = [1,2;2,3;3,4;4,1;
               5,6;6,7;7,8;8,5;
               2,3;3,7;7,6;6,2;
               1,4;4,8;8,5;5,1;
               1,2;2,6;6,5;5,1;
               3,4;4,8;8,7;7,3];
         
end

function d = sdfCube(P,x1,x2,y1,y2,z1,z2)
d = [x1-P(:,1) + 1e-12, P(:,1)-x2, y1-P(:,2)+1e-12, P(:,2)-y2, z1-P(:,3)+1e-12, P(:,3)-z2];
d = [d,max(d,[],2)];
end
