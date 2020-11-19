function Dist = imsdf(im,p)
D1 = bwdist(im);
D2 = bwdist(-1.0*(im-255));
% 
D = D1 - D2;
[X,Y] = meshgrid(1:size(im,2),1:size(im,1));

V = double(interp2(X,Y,D,p(:,1),p(:,2)));
V(isnan(V(:))) = 100;
% 
Dist = [V,V,V,V];

% B = bwboundaries(im);
% X = vertcat(B{:});
% 
% P(:,2) = X(:,1) - 1;
% P(:,1) = X(:,2) - 1;
% 
% [~, d] = knnsearch(P,Y);
% 
% [x,y] = meshgrid(1:size(im,2),1:size(im,1));
% Sign = double(interp2(x,y,(-2*(im==255)+1),Y(:,1),Y(:,2)));
% 
% Dist = d;
end

