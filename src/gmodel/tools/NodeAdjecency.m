%------------------------------------- GENERATE ADJECENCY MATRIX FROM FACES
function A = NodeAdjecency(face)
n = max(cellfun(@max,face));      
A = sparse(n,n);
null = sparse(n,n);

for ii = 1:length(face)
    poly = uint16(face{ii});
    [i,j] = PermutationSet(numel(poly));
    B = null;
    B(poly(i),poly(j)) = 1;
    A = A + B;
end

A = double(A>0);
end

%---------------------------------------------------------- PERMUTATION SET
function [i,j] = PermutationSet(n)
i = transpose(kron(1:n,ones(1,n)));
j = transpose(kron(ones(1,n),1:n));
end