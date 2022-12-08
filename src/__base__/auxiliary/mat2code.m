function mat2code(Matrix,Name,Size)
% MAT2CODE.m generate matlab code for symbolic matrix.

for i = 1:Size(1)
    for j = 1:Size(2)
        disp([Name, '(', num2str(i), ',', num2str(j), ') = ', char(Matrix(i, j)), ';']);
    end
end

end