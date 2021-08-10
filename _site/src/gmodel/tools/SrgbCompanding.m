%------------------------------------------- CONVERT/INVERT sRGB COMPANDING
function ColorNew = SrgbCompanding(Color,Inverse)
r = Color(1);
g = Color(2);
b = Color(3);

if Inverse == 1
    if (r > 0.04045), r = ((r+0.055)/1.055)^2.4; else, r = r / 12.92; end
    if (g > 0.04045), g = ((g+0.055)/1.055)^2.4; else, g = g / 12.92; end
    if (b > 0.04045), b = ((b+0.055)/1.055)^2.4; else, b = b / 12.92; end
else
    if (r > 0.0031308),r = 1.055*(r)^(1/2.4)-0.055; else, r = r * 12.92;end
    if (g > 0.0031308),g = 1.055*(g)^(1/2.4)-0.055; else, g = g * 12.92;end
    if (b > 0.0031308),b = 1.055*(b)^(1/2.4)-0.055; else, b = b * 12.92;end
end

ColorNew = [r, g, b];
end