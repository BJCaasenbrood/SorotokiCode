function rgbImage = colorshift(rgbImage,alpha)

%rgbImage = rgb2lin(rgbImage);
hsvImage = rgb2hsv(rgbImage);
hsvImage(:,:,1) = rem(hsvImage(:,:,1) * alpha, 1);
rgbImage = hsv2rgb(hsvImage);
%rgbImage = lin2rgb(rgbImage);
rgbImage = clamp(rgbImage*255,0,255);

end

