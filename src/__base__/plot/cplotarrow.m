function h = cplotarrow(x,y,c,map,varargin)
% orginal code by Michael Heidingsfeld 
% Plot data:
if ~(all(size(x) == size(y)) && all(size(x) == size(c)))
    error('Vectors X,Y and C must be the same size');
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
    h(k) = plot(X,Y,'Color',map(k,:),varargin{:});
    
    if mod(k,6) == 0 
        X_ = X(~isnan(X)); nX = diff(X_);
        Y_ = Y(~isnan(Y)); nY = diff(Y_);
        nrm = [nX(end),nY(end)];
        nrm = nrm./norm(nrm);
        %NRM = [NRM;nrm];
        V   = [X_(end),Y_(end)];
        
        [V,F] = arrowhead(V,atan2(nrm(1),nrm(2)));
        
        patch('Vertices',V,'Faces',F,'EdgeColor','None',...
            'FaceColor',map(k,:));
    end
    hold on;
end
if status == 0, hold off; end   % reset hold status
end

function [V,F] = arrowhead(X,th)
a = 0.22;
b = 0.77;
V = 0.175*[-a,-0.2*b;
     0,0;
     a,-0.2*b;
     0,b]; 
 
th = th - 0.1;
R = [cos(th),sin(th);-sin(th),cos(th)];
V = (R*V.').';
 
V(:,1) = V(:,1) + X(1);
V(:,2) = V(:,2) + X(2);
 
F = [1,2,3,4];
end

