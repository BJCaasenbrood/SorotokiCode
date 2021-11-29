function fem = TabulateShapeFunctionsC3T10(fem)

ElemNNode    = cellfun(@length,fem.Element); 
fem.ShapeFnc = cell(max(ElemNNode),1);

[W,Q] = Tetrahedron();

fem.ShapeFnc{10}.W = W;
fem.ShapeFnc{10}.N = zeros(4,1,size(W,1));
fem.ShapeFnc{10}.dNdxi = zeros(4,3,size(W,1));
fem.ShapeFnc{10}.Q = Q;

for q = 1:size(W,1)
    [N,dNdxi] = TetraShapeFnc(Q(q,:));
    fem.ShapeFnc{4}.N(:,:,q) = N;
    fem.ShapeFnc{4}.dNdxi(:,:,q) = dNdxi;
end

end

function [weight,point] = Tetrahedron

xa = [0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105]; 
ya = [0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685];
za = [0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105];
wt = [0.2500000000000000, 0.2500000000000000, 0.2500000000000000, 0.2500000000000000]/6;

xw=[xa(1),ya(1),za(1),wt(1);
    xa(2),ya(2),za(2),wt(2);
    xa(3),ya(3),za(3),wt(3);
    xa(4),ya(4),za(4),wt(4)];

point  = xw(:,1:3);
weight = xw(:,4);
end

function [N, dNdxi] = TetraShapeFnc(s)
 
N     = zeros(10,1);
dNdxi = zeros(10,3);
N(1) = 1 - s(1) - s(2) - s(3);
N(2) = s(1);
N(3) = s(2);
N(4) = s(3);

dNdxi(1,1) = -1;
dNdxi(2,1) = 1;
dNdxi(3,1) = 0;
dNdxi(4,1) = 0;
dNdxi(1,2) = -1;
dNdxi(2,2) = 0;
dNdxi(3,2) = 1;
dNdxi(4,2) = 0;
dNdxi(1,3) = -1;
dNdxi(2,3) = 0;
dNdxi(3,3) = 0;
dNdxi(4,3) = 1;

end

