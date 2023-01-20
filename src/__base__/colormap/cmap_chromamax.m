function colmap = chromamax(alpha)
    I = imread('chromamax.jpeg');
    id = round(500*clamp(alpha,0,1));
    colmap = squeeze(I(id,:,:));
end

