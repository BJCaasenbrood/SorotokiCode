function h = cplot(x,y,c,map,varargin)
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
    hold on;
end

if status == 0
    hold off; % reset hold status
end   

end

