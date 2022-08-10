function colmap = cmap_goldgreen(varargin)

if isempty(varargin)
    n = 200;
else
    n = varargin{1};
end

colmap = [color_gold;
         color_gold_light
         color_gold_lighter;
         color_gold_lightest;
         color_green_lightest
         color_green_lighter;
         color_green_light;
         color_green];

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