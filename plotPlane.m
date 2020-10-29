function plotPlane(pts)
[~,V,p] = affine_fit(pts);

[S1,S2] = meshgrid([-0.5 0 0.5]);

%generate the pont coordinates
X = p(1)+([S1(:) S2(:)]*V(1,:)');
Y = p(2)+([S1(:) S2(:)]*V(2,:)');
Z = p(3)+([S1(:) S2(:)]*V(3,:)');

[I,~] = imread('checker.jpg');

%plot the plane
warpim(reshape(X,3,3),reshape(Y,3,3),reshape(Z,3,3),I);
hold on;

patch('Faces',[1,2,3,4],'Vertices',pts,...
    'Linewidth',1.5,'linestyle','-','FaceColor','none',...
    'EdgeColor',[1 1 1]*0.5);

function [n,V,p] = affine_fit(X)
p = mean(X,1);
R = bsxfun(@minus,X,p);
[V,~] = eig(R'*R);
n = V(:,1);
V = V(:,2:end);
end


end

