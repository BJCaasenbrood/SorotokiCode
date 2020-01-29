function T = TensorOperation(A,B,Arg)
% tic
% [id, set, W] = Tensor4IdSymmetric;

%T = zeros(length(id)^2,1);

if strcmp(Arg,'x') % kronecker product
    T = TensorProduct(A,B);
%     for kk = 1:length(set)
%         row = set{kk}(1);
%         col = set{kk}(2);
%         i = id(row,1);
%         j = id(row,2);
%         k = id(col,1);
%         l = id(col,2);
%         Aij = A(i,j);
%         Bkl = B(k,l);
%         T(kk) = Aij*Bkl;
%     end
elseif strcmp(Arg,'xt') % symmetric kronecker product
    T = SymmetricTensorProduct(A,B);
%     for kk = 1:length(set)
%         row = set{kk}(1);
%         col = set{kk}(2);
%         i = id(row,1);
%         j = id(row,2);
%         k = id(col,1);
%         l = id(col,2);
%         Aik = A(i,k);
%         Ail = A(i,l);
%         Bjk = B(j,k);
%         Bjl = B(j,l);
%         
%         T(kk) = 0.5*(Aik*Bjl + Ail*Bjk);
%     end
end
% 
% T = reshape(T,length(id),length(id));
% T = T.*W;
% toc
end

function T = TensorProduct(A,B)

a11 = A(1,1); a12 = A(1,2); a13 = A(1,3);  
a21 = A(2,1); a22 = A(2,2); a23 = A(2,3);  
a31 = A(3,1); a32 = A(3,2); a33 = A(3,3);  

b11 = B(1,1); b12 = B(1,2); b13 = B(1,3);  
b21 = B(2,1); b22 = B(2,2); b23 = B(2,3);  
b31 = B(3,1); b32 = B(3,2); b33 = B(3,3);  

T = [   a11*b11,   a22*b11,   a33*b11, 2*a12*b11, 2*a23*b11, 2*a13*b11;
   a11*b22,   a22*b22,   a33*b22, 2*a12*b22, 2*a23*b22, 2*a13*b22;
   a11*b33,   a22*b33,   a33*b33, 2*a12*b33, 2*a23*b33, 2*a13*b33;
 2*a11*b12, 2*a22*b12, 2*a33*b12, 4*a12*b12, 4*a23*b12, 4*a13*b12;
 2*a11*b23, 2*a22*b23, 2*a33*b23, 4*a12*b23, 4*a23*b23, 4*a13*b23;
 2*a11*b13, 2*a22*b13, 2*a33*b13, 4*a12*b13, 4*a23*b13, 4*a13*b13];
end

function T = SymmetricTensorProduct(A,B)

a11 = A(1,1); a12 = A(1,2); a13 = A(1,3);  
a21 = A(2,1); a22 = A(2,2); a23 = A(2,3);  
a31 = A(3,1); a32 = A(3,2); a33 = A(3,3);  

b11 = B(1,1); b12 = B(1,2); b13 = B(1,3);  
b21 = B(2,1); b22 = B(2,2); b23 = B(2,3);  
b31 = B(3,1); b32 = B(3,2); b33 = B(3,3);  

T = [      a11*b11,           a21*b21,           a31*b31,             2*a11*b21,             2*a21*b31,             2*a11*b31;
           a12*b12,           a22*b22,           a32*b32,             2*a12*b22,             2*a22*b32,             2*a12*b32;
           a13*b13,           a23*b23,           a33*b33,             2*a13*b23,             2*a23*b33,             2*a13*b33;
 a11*b12 + a12*b11, a21*b22 + a22*b21, a31*b32 + a32*b31, 2*a11*b22 + 2*a12*b21, 2*a21*b32 + 2*a22*b31, 2*a11*b32 + 2*a12*b31;
 a12*b13 + a13*b12, a22*b23 + a23*b22, a32*b33 + a33*b32, 2*a12*b23 + 2*a13*b22, 2*a22*b33 + 2*a23*b32, 2*a12*b33 + 2*a13*b32;
 a11*b13 + a13*b11, a21*b23 + a23*b21, a31*b33 + a33*b31, 2*a11*b23 + 2*a13*b21, 2*a21*b33 + 2*a23*b31, 2*a11*b33 + 2*a13*b31];
end