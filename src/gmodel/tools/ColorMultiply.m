%------------------------------------------------------ MULTIPLY RBG COLORS
function N = ColorMultiply(N1,N2)
N = zeros(size(N1,1),3);

for i = 1:length(N)
    N(1) = N1(1)*N2(1);
    N(2) = N1(2)*N2(2);
    N(3) = N1(3)*N2(3);
end
end