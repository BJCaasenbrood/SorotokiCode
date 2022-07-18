% Y = ANNIHIL(X) computes the annihilator of X such that Y*X = 0.
function Y = annihil(X)
Y0 = zeros(size(X));
Y  = null(X.').';

% if rank(Y) == 0
%     Y = Y0;
% end
end

