function Y = invchol(X)
% Matrix Inversion using Cholesky Decomposition
%
% Finds the inverse of the matrix X, given its (lower triangular) Cholesky
% Decomposition; i.e. X = LL', according to the paper 'Matrix Inversion
% Using Cholesky Decomposition', Aravindh Krishnamoorthy, Deepak Menon,
% arXiv:1111.4144.
%

% Version 0.1, 2013-05-25, Aravindh Krishnamoorthy
% e-mail: aravindh.k@ieee.org
L = chol(X,'lower');
N = size(L, 1);
U = L\eye(N);

Y = L'\U;
% Y = zeros(N, N) ;
% R = L' ;
% S = (diag(diag(R)))\eye(N);
% 
% for j=N:-1:1
%     for i=j:-1:1
%         Y(i,j) = S(i,j) - R(i,i+1:end)*Y(i+1:end,j) ;
%         Y(i,j) = Y(i,j)/R(i,i) ;
%         Y(j,i) = conj(Y(i,j)) ;
%     end
% end
end
