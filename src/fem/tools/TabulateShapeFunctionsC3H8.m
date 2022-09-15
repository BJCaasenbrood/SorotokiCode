function fem = TabulateShapeFunctionsC3H8(fem)

ElemNNode    = cellfun(@length,fem.Element); 
fem.ShapeFnc = cell(max(ElemNNode),1);

[W,Q] = Hexahedron();

fem.ShapeFnc{8}.W = W;
fem.ShapeFnc{8}.N = zeros(8,1,size(W,1));
fem.ShapeFnc{8}.dNdxi = zeros(8,3,size(W,1));
fem.ShapeFnc{8}.Q = Q;
fem.ShapeFnc{8}.fnc = @(x) HexaShapeFnc(x);
fem.ShapeFnc{8}.Xi = [-1  1  1 -1 -1  1  1 -1;
                      -1 -1  1  1 -1 -1  1  1;
                      -1 -1 -1 -1  1  1  1  1].';

for q = 1:size(W,1)
    [N,dNdxi] = HexaShapeFnc(Q(q,:));
    fem.ShapeFnc{8}.N(:,:,q) = N;
    fem.ShapeFnc{8}.dNdxi(:,:,q) = dNdxi;
end

end

function [weight,point] = Hexahedron

xw=...
[-0.5773502691896  -0.5773502691896  -0.5773502691896  1.00000000000000
 -0.5773502691896  -0.5773502691896   0.5773502691896  1.00000000000000
 -0.5773502691896   0.5773502691896  -0.5773502691896  1.00000000000000
 -0.5773502691896   0.5773502691896   0.5773502691896  1.00000000000000
  0.5773502691896  -0.5773502691896  -0.5773502691896  1.00000000000000
  0.5773502691896  -0.5773502691896   0.5773502691896  1.00000000000000
  0.5773502691896   0.5773502691896  -0.5773502691896  1.00000000000000
  0.5773502691896   0.5773502691896   0.5773502691896  1.00000000000000];

point  = xw(:,1:3);
weight = xw(:,4);
end

function [N, dNdxi] = HexaShapeFnc(s)

Xez = [-1  1  1 -1 -1  1  1 -1;
       -1 -1  1  1 -1 -1  1  1;
       -1 -1 -1 -1  1  1  1  1];
  
N     = zeros(8,1);
dNdxi = zeros(8,3);

for ii=1:8
    Xp = Xez(1,ii);
    Yp = Xez(2,ii);
    Zp = Xez(3,ii);
    
    Xi0 = [1+s(1)*Xp 1+s(2)*Yp 1+s(3)*Zp];
    
    N(ii) = 0.125*Xi0(1)*Xi0(2)*Xi0(3);
    dNdxi(ii,1) = 0.125*Xp*Xi0(2)*Xi0(3);
    dNdxi(ii,2) = 0.125*Yp*Xi0(1)*Xi0(3);
    dNdxi(ii,3) = 0.125*Zp*Xi0(1)*Xi0(2);
end

end

