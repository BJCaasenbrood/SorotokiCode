%-------------------------------------------- CONVERT NORMALS TO RGB COLORS
function [N, Material] = Normal2RGB(Domain,N)    

if isempty(Domain('Material')), Material = 'Porcelain';
else, Material = Domain('Material');
end

if size(N,2) < 3, N = repmat(N,1,3); end

C1 = ColorScheme(Material,1)*0.95;
C2 = ColorScheme(Material,2)*0.95;
C3 = ColorScheme(Material,3)*0.95;

Ntmp = zeros(length(N),3);
if strcmp(Material,'Bump')
    N(:,1) = 0.5 + 0.5*(N(:,1));
    N(:,2) = 0.5 + 0.5*(N(:,2));
    N(:,3) = 0.5 + 0.5*softabs(N(:,3));
else
N(:,1) = 0.1 + 0.6*softabs(N(:,1));
N(:,2) = 0.1 + 0.6*softabs(N(:,2));
N(:,3) = 0.3 + 0.6*(N(:,3));
end

CB = ColorMultiply(ColorMultiply(C1,C2),C3);
for ii = 1:length(N)
    tmp = AffineColorMix(C1,C2,C3,N(ii,:));
    Ntmp(ii,:) = [max(CB(1),tmp(1)),...
        max(CB(2),tmp(2)),...
        max(CB(3),tmp(3))];
end
N = Ntmp;
end