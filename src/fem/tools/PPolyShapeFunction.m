%------------------------------------------------ POLYGONAL SHAPE FUNCTIONS
function [N,dNdxi] = PolyShapeFunction(nn,xi)
N = zeros(nn,1);
alpha  = zeros(nn,1);
dNdxi  = zeros(nn,2);
dalpha = zeros(nn,2);
sum_alpha = 0;
sum_dalpha = zeros(1,2);
A  = zeros(nn,1);
dA = zeros(nn,2);
[p,Tri] = PolyTrnglt(nn,xi);

for i=1:nn
    sctr = Tri(i,:);
    pT = p(sctr,:);
    A(i,1) = 1/2*det([pT,ones(3,1)]);
    dA(i,1) = 1/2*(pT(3,2)-pT(2,2));
    dA(i,2) = 1/2*(pT(2,1)-pT(3,1));
end

A  = [A(nn,:);A];
dA = [dA(nn,:);dA];

for i=1:nn
    alpha(i,1) = 1/(A(i)*A(i+1));
    dalpha(i,1) = -alpha(i)*(dA(i,1)/A(i)+dA(i+1,1)/A(i+1));
    dalpha(i,2) = -alpha(i)*(dA(i,2)/A(i)+dA(i+1,2)/A(i+1));
    sum_alpha = sum_alpha + alpha(i);
    if i == 1, sum_dalpha(1:2) = dalpha(i,1:2);
    else, sum_dalpha(1:2) = sum_dalpha(1:2) + dalpha(i,1:2);
    end
end

for i=1:nn
    N(i) = alpha(i)/sum_alpha;
    dNdxi(i,1:2) = (dalpha(i,1:2)-N(i)*sum_dalpha(1:2))/sum_alpha;
end

end

%---------------------------------------------------- POLYGON TRIANGULATION
function [p,Tri] = PolyTrnglt(nn,xi)
p = [cos(2*pi*((1:nn))/nn);
     sin(2*pi*((1:nn))/nn)]';

p = [p; xi];

Tri = zeros(nn,3);
Tri(1:nn,1) = nn+1;
Tri(1:nn,2) = 1:nn;
Tri(1:nn,3) = 2:nn+1;
Tri(nn,3) = 1;
end
