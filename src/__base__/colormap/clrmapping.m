function colmap = clrmapping(colmap,arg)

if arg > 1
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,arg);
    colmap = interp1(x,colmap,xq);
elseif arg < 0
    if abs(arg) == 1, arg = 1e3; end
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,abs(arg));
    colmap = interp1(x,colmap,xq);
    colmap = flipud(colmap);
elseif arg == 0
    A = flipud(colmap(1:50,:));
    B = flipud(colmap(51:end,:));
    colmap = [A;flipud(A);flipud(B);B];
    %colmap = [flipud(A);flipud(B)];
end

end