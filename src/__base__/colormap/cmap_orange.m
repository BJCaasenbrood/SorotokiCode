function colmap = cmap_orange(varargin)

if isempty(varargin)
    C = color_orange;
    n = 100;
else
    if numel(varargin) == 1 && numel(varargin{1}) == 1
        if numel(varargin{1}) == 1
            n = varargin{1};
            C = col(1);
        elseif numel(varargin{1}) == 3
            C = varargin{1};
            n = 100;
        end
        
    elseif numel(varargin) == 2
        if numel(varargin{1}) == 1 && numel(varargin{2}) == 3
            n = varargin{1};
            C = varargin{2};
        elseif numel(varargin{2}) == 1 && numel(varargin{1}) == 3
            C = varargin{1};
            n = varargin{2};
        end
    end
end

WW = [0.8929    0.8929    0.8929];


if nargin>0
    colmap = linspacen(WW.',C.',100).';
    colmap = clrmapping(colmap,n);
else
    colmap = linspacen(WW.',C.',n).';
end

end

function colmap = clrmapping(colmap,arg)

if arg > 1
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,arg);
    colmap = interp1(x,colmap,xq);
elseif arg < 0
    if abs(arg) == 1, arg = 100; end
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,abs(arg));
    colmap = interp1(x,colmap,xq);
    colmap = flipud(colmap);
elseif arg == 0
    colmap = [colmap;flipud(colmap)];
end

end