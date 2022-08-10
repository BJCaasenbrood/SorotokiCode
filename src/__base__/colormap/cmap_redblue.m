function colmap = cmap_redblue(varargin)

if isempty(varargin)
    n = 200;
else
    n = varargin{1};
end

colmap = [color_primary_alt_darkest;
         color_primary_alt_dark;
         color_primary_alt_light;
         color_primary_alt_lightest;
         color_secondary_lightest
         color_secondary_light;
         color_secondary_dark;
         color_secondary_darkest];

colmap = clrmapping(colmap,n);

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