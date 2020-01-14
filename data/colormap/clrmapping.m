function colmap = clrmapping(colmap,arg)

if arg > 1
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,n);
    colmap = interp1(x,colmap,xq);
elseif arg < 0
    colmap = flipud(colmap);
elseif arg == 0
    colmap = [colmap;flipud(colmap)];
end

end