%------------------------------------------------ AFFINE COLOR COMBINATIONS
function Color = AffineColorMix(C1,C2,C3,Mix)
a1 = Mix(1);
a2 = Mix(2);
a3 = Mix(3);

Color1 = SrgbCompanding(C1,1);
Color2 = SrgbCompanding(C2,1);
Color3 = SrgbCompanding(C3,1);

r = Color1(1)*a1 + Color2(1)*a2 + Color3(1)*a3;
g = Color1(2)*a1 + Color2(2)*a2 + Color3(2)*a3;
b = Color1(3)*a1 + Color2(3)*a2 + Color3(3)*a3;

Color = [r,g,b];

Color = SrgbCompanding(Color,0);
end