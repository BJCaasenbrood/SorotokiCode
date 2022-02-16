function colmap = diffmap(h,n)

HSV = 1.05*rgb2hsv(noir);
HSV(:,1) = h;

colmap = hsv2rgb(HSV);

if nargin>0
    colmap = clrmapping(colmap,n);
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