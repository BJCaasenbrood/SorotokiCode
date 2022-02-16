function mcap = diffuse(varargin)

if isempty(varargin), h = 0.5;
else, h = varargin{1};
end

I = 1.15*rgb2hsv(imread('classicblue.jpg'));
I(:,:,1) = ones(size(I,1),size(I,2))*h;
%I(:,:,2) = I(:,:,2);

mcap = uint8(round(hsv2rgb(I)*255));
end


