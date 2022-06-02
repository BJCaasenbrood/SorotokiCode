function h = cplot3(x,y,z,c,map,varargin)
% orginal code by Michael Heidingsfeld 
% Plot data:
if ~(all(size(x) == size(y)) && all(size(x) == size(c)))
    error('Vectors X,Y,Z and C must be the same size');
end
N = size(map,1);
cmax = max(c);
cmin = min(c);
cint = (cmax-cmin)/N;
indices = 1:length(c);
status = ishold;                % save hold status
%NRM = [];
%V   = [];
for k = 1:N
    ii = logical(c >= cmin+k*cint) + logical(c <= cmin+(k-1)*cint);
    jj = ones(size(ii)); jj(1:end-1) = ii(2:end);
    ii = logical(ii .* jj);
    X = x; X(indices(ii)) = NaN;
    Y = y; Y(indices(ii)) = NaN;
    Z = z; Z(indices(ii)) = NaN;
    h(k) = plot3(X,Y,Z,'Color',map(k,:),varargin{:});
    hold on;
%     if mod(k,6) == 0 
%         X_ = X(~isnan(X)); nX = diff(X_);
%         Y_ = Y(~isnan(Y)); nY = diff(Y_);
%         Z_ = Z(~isnan(Z)); nZ = diff(Z_);
%         nrm = [nX(end),nY(end),nZ(end)];
%         nrm = nrm./norm(nrm);
%         %NRM = [NRM;nrm];
%         V   = [X_(end),Y_(end),Z_(end)];
%         
%         R = vrrotvec2mat(vrrotvec([0,1,0].',nrm(:)));
%         
%         [V,F] = arrowhead(V,R);
%         
%         patch('Vertices',V,'Faces',F,'EdgeColor','None',...
%             'FaceColor',map(k,:));
%     end
%     hold on;
end
if status == 0, hold off; end   % reset hold status
end

function [V,F] = arrowhead(X,R)
a = 0.22;
b = 0.77;
V = 0.00375*[
    -a,-0.2*b,0;
     0,0,0;
     a,-0.2*b,0;
     0,b,0
    0,-0.2*b,-a;
     0,0,0;
    0,-0.2*b,a;
     0,b,0]; 
 
V = (R*V.').';
 
V(:,1) = V(:,1) + X(1);
V(:,2) = V(:,2) + X(2);
V(:,3) = V(:,3) + X(3);
 
F = [1,2,3,4;5,6,7,8];
end



