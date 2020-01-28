function Img = ImageColormix(img,Color)
Img = zeros(size(img));

for ii = 1:size(img,1)
    for jj = 1:size(img,2)
        
       C1 = reshape(double(img(ii,jj,:)),[3,1]);
       tmp = AffineMixing(C1,Color,[0.5,0.5]);
       Img(ii,jj,:) = round(tmp);
    end
end


end

%------------------------------------------------ AFFINE COLOR COMBINATIONS
function Color = AffineMixing(C1,C2,Mix)
a1 = Mix(1);
a2 = Mix(2);

Color1 = SrgbCompanding(C1,0);
Color2 = SrgbCompanding(C2,0);

r = Color1(1)*a1 + Color2(1)*a2;
g = Color1(2)*a1 + Color2(2)*a2;
b = Color1(3)*a1 + Color2(3)*a2;

Color = [r,g,b];

Color = SrgbCompanding(Color,1);
end