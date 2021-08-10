function Img = ImageColormix(img,Color)
Img = zeros(size(img));

if max(abs(Color)) <= 1, Color = Color*255; end

for ii = 1:size(img,1)
    for jj = 1:size(img,2)
       C1 = reshape(double(img(ii,jj,:)),[3,1]);
       tmp = AffineMixing(C1,Color,[0.3,0.7]);
       Img(ii,jj,:) = round(tmp);
    end
end

end

%------------------------------------------------ AFFINE COLOR COMBINATIONS
function Color = AffineMixing(C1,C2,Mix)
a = Mix(1);
b = Mix(2);

Color1 = SrgbCompanding(C1,0);
Color2 = SrgbCompanding(C2,0);

r = Color1(1)*a + Color2(1)*b;
g = Color1(2)*a + Color2(2)*b;
b = Color1(3)*a + Color2(3)*b;
% 
% r = Color1(1)*Color2(1)/(255);
% g = Color1(2)*Color2(2)/(255);
% b = Color1(3)*Color2(3)/(255);
% 
Color = [r,g,b];
% 
% Color = 255*mix(Color1/255,Color2/255);

Color = SrgbCompanding(Color,1);
end

function c = mix(a,b)
c = a.*b;

end