%------------------------------------------ CONVERT RGB COLORS TO GRAYSCALE
function N = RGB2Gray(N)

Clin = 0.2126*N(:,1) + 0.7152*N(:,2) + 0.0722*N(:,3);
Csrgb = zeros(size(N,1),1);
for i = 1:size(N,1)
if Clin(i) <= 0.031308
    Csrgb(i,:) = 12.92*Clin(i);
else
    Csrgb(i,:) = 1.055*Clin(i).^(1/2.4) - 0.055;
end
end

N(:,1) = Csrgb*0.98;
N(:,3) = Csrgb*0.98;
N(:,2) = Csrgb*0.98;
end