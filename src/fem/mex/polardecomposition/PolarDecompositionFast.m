%#codegen
function [R,S,V] = PolarDecompositionFast(F)

C = F'*F;
[Q0, lambdasquare] = eig(C);

lambda = sqrt(diag((lambdasquare))); 
Uinv   = repmat(1./lambda',size(F,1),1).*Q0*Q0';

R = F*Uinv;
S = R'*F;
V = F*R';

end