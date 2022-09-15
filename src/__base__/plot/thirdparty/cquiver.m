function quiv = cquiver(varargin)
% myquiver quiver with color
%
% <SYNTAX>
%   cquiver(varargin);
%
% <DESCRIPTION>
% cquiver(x, y, u, v) plots quiver for at given points (x, y)
% with the velocity (u, v)
%
% quiv = cquiver(x, y, u, v) returns quiver structure with two properties,
% Head and Tail.  Here, Head and Tail are the patch.
%
% quiv = cquiver(x, y, u, v, 'cdata', 'angle') set color based on angle.
%
% quiv = cquiver(x, y, u, v, 'sampling', h) plots quiver plot
% at uniform grid where grid size is h.
% Default sampling method is average.
% Sampling method can be changed to 'fastest' by using
% quiv = cquiver(x, y, u, v, 'sampling', h, 'samplingmethod', 'fastest')
% In this case, the fastest one will be chosen from the rectangle
% (X(i)-h/2, X(i)+h/2) x (Y(i)-h/2, Y(i)+h/2).
%
% See also, quiver
% Copyright 2019 Dohyun Kim / CC BY-NC
% Contact: kim92n@gmail.com
% Developed using MATLAB.ver 9.6 (R2019a) on Microsoft Windows 10 Enterprise
%%
if ~isa(varargin{1}, 'matlab.graphics.axis.Axes')
    cquiver(gca, varargin{:});
    return;
end
ax = varargin{1};
x = varargin{2};
y = varargin{3};
u = varargin{4};
v = varargin{5};
names = varargin(6:2:end);
vals  = varargin(7:2:end);
names = cellfun(@lower, names, 'un', 0);
% set arrow head scale
isHeadScaleFactor = strcmpi(names, 'headscale');
if any(isHeadScaleFactor)
    headScaleFactor = vals(isHeadScaleFactor);
    headScaleFactor = headScaleFactor{end};
else
    headScaleFactor = 0.3;
end
names(isHeadScaleFactor) = [];
vals(isHeadScaleFactor) = [];
% get arrow length
lims = axis(ax);
isScale = strcmpi(names, 'scale');
if any(isScale)
    scale = vals(isScale);
    scale = scale{end};
else
    scale = gmean(diff(reshape(lims,2,[])))./sqrt(numel(x))*2;
end
names(isScale) = [];
vals(isScale) = [];
% set arrow head scale
isCData = strcmpi(names, 'cdata');
if any(isCData)
    CData = vals(isCData);
    CData = CData{end};
    switch CData
        case 'angle'
            isAngular = true;
        case 'magnitude'
            isAngular = false;
        otherwise
            error('CData should be ''angle'' or ''magnitude''')
    end
else
    isAngular = false;
end
names(isCData) = [];
vals(isCData) = [];
% set uniform sampling
isSampling = strcmpi(names, 'sampling');
if any(isSampling)
    h = vals{isSampling(end)};
    isSamplingMethod = strcmpi(names, 'samplingmethod');
    names(isSampling) = [];
    vals(isSampling) = [];
    xlim = [min(x(:)), max(x(:))];
    ylim = [min(y(:)), max(y(:))];
    if ~isempty(xlim) && ~isempty(ylim)
        [X,Y] = ndgrid(linspace(xlim(1), xlim(2), round(diff(xlim)/h)), linspace(ylim(1), ylim(2), round(diff(ylim)/h)));
        X = X(:); Y = Y(:);
        if any(isSamplingMethod)
            SamplingMethod = vals(isSamplingMethod);
            SamplingMethod = lower(SamplingMethod{end});
        else
            SamplingMethod = 'average';
        end
        xc = zeros(size(X)); yc = zeros(size(X));
        uc = zeros(size(X)); vc = zeros(size(Y));
        hasPoint = true(size(X));
        switch SamplingMethod
            case 'fastest'
                vel = sqrt(u.^2 + v.^2);
                for i = 1 : numel(X)
                    idx = max(abs(x-X(i)),abs(y-Y(i))) < h/2;
                    if nnz(idx)
                        [~,ii] = max(vel(idx));
                        idx(idx) = accumarray(ii,true,[nnz(idx),1]);
                        uc(i) = u1(idx); vc(i) = u2(idx);
                        xc(i) = x(idx); yc(i) = Y(idx);
                    else
                        hasPoint(i) = false;
                    end
                end
            case 'average'
                for i = 1 : numel(X)
                    idx = max(abs(x-X(i)), abs(y-Y(i))) < h/2;
                    if nnz(idx)
                        xc(i) = mean(x(idx)); yc(i) = mean(y(idx));
                        uc(i) = mean(u(idx)); vc(i) = mean(v(idx));
                    else
                        hasPoint(i) = false;
                    end
                end
        end
        x = xc(hasPoint); y = yc(hasPoint);
        u = uc(hasPoint); v = vc(hasPoint);
        names(isSamplingMethod) = [];
        vals(isSamplingMethod) = [];
    end
end
%% QUIVER
vel = sqrt(u.^2 + v.^2);
maxMag = max(vel(:));
u = u./maxMag;
v = v./maxMag;
% compute arrow end point
x1 = x + u*scale;
y1 = y + v*scale;
headscale = scale*headScaleFactor;
% get normal vector
nx =  v;
ny = -u;
% get arrow mid point
mx = x + u*(scale-headscale);
my = y + v*(scale-headscale);
xl = mx + nx*headscale/3;
xr = mx - nx*headscale/3;
yl = my + ny*headscale/3;
yr = my - ny*headscale/3;
args = [names;vals];
if isAngular
    t = cart2pol(u(:), v(:));
    quiv.Head = patch('XData', [x1(:).';xl(:).';xr(:).'], 'YData', [y1(:).';yl(:).';yr(:).'], 'CData', repmat(t(:).',3,1), 'FaceColor' ,'flat', 'EdgeColor', 'none', args{:});
    isFaceOpt = contains(names, 'face');
    args(:,isFaceOpt) = [];
    quiv.Tail = patch('XData', [x(:).';x1(:).'], 'YData', [y(:).';y1(:).'], 'CData', repmat(t(:).',2,1), 'FaceColor' ,'none', 'EdgeColor', 'flat', args{:});
else
    quiv.Head = patch('XData', [x1(:).';xl(:).';xr(:).'], 'YData', [y1(:).';yl(:).';yr(:).'], 'CData', repmat(vel(:).',3,1), 'FaceColor' ,'flat', 'EdgeColor', 'none', args{:});
    isFaceOpt = contains(names, 'face');
    args(:,isFaceOpt) = [];
    quiv.Tail = patch('XData', [x(:).';x1(:).'], 'YData', [y(:).';y1(:).'], 'CData', repmat(vel(:).',2,1), 'FaceColor' ,'none', 'EdgeColor', 'flat', args{:});
end
