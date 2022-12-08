% Y = ANNIHIL(X) computes the annihilator of X such that Y*X = 0.
function Y = annihil(X)
Y  = null(X.').';
end

