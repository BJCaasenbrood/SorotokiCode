% D1 = bwdist(I);
% D2 = bwdist(-1*(I-255));
% 
% D = D1 - D2;
% [X,Y] = meshgrid(1:300,1:600);
% 
% figure(1);
% surf(X,Y,D,'linestyle','none');
% 
% Vq = interp2(X,Y,D,300,606)

% id = find(Vq < 0);
% 
% figure(2);
% plot(X(id),Y(id),'.');

% Dist = @(p) imsdf(I,p);

msh = Mesh(I,'BdBox',[0,size(I,2),0,size(I,1)],'NElem',3e3);
msh = msh.generate();