function mcap = diffuse(h)
I = 1.15*rgb2hsv(imread('classicblue.jpg'));
I(:,:,1) = ones(size(I,1),size(I,2))*h;
I(:,:,2) = I(:,:,2);

mcap = round(hsv2rgb(I)*255);
end


