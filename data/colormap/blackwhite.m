function colmap = blackwhite(n)

if nargin>0
    if n < 1, colmap = gray(200);
    else, colmap = gray(n); end
    colmap = clrmapping(colmap,n);
else
    colmap = gray(200);
end

end

