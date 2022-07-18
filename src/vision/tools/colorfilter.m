function [I,mask] = colorfilter(image, range)
    % RGB to HSV conversion
    I = rgb2hsv(image);         
    
    % Normalization range between 0 and 1
    range = range./360;
    
    % Mask creation
    if(size(range,1) > 1), error('Error. Range matriz has too many rows.'); end
    if(size(range,2) > 2), error('Error. Range matriz has too many columns.'); end
    if(range(1) > range(2))
        % Red hue case
        mask = (I(:,:,1)>range(1) & (I(:,:,1)<=1)) + (I(:,:,1)<range(2) & (I(:,:,1)>=0));
    else
        % Regular case
        mask = (I(:,:,1)>range(1)) & (I(:,:,1)<range(2));
    end
    
    % Saturation is modified according to the mask
    I(:,:,2) = mask .* I(:,:,2);
    
    % HSV to RGB conversion
    I = hsv2rgb(I);
    
end