function mat2code(Matrix,Name,Size)
%diary([Name,'.txt'])
for i = 1:Size(1)
for j = 1:Size(2)
disp([Name,'(',num2str(i),',',num2str(j),') = ',char(Matrix(i,j)),';']);
end
end
%diary off
end