function A = overwrite_colors(colmap)
A = [];

if nargin < 1
    for ii = 1:12
        A = vappend(A,col(ii)); 
    end
else
    A = colmap; 
end

set(gca, 'ColorOrder', A, 'NextPlot', 'replacechildren');
end

