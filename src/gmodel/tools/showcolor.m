function showcolor(varargin)
for ii = 1:nargin
    dx = ii-1;
    ds = 0.05;
    v = [dx,0;dx+1-ds,0;dx+1-ds,1;dx,1];
    f = [1,2,3,4];
    patch('Vertices',v,'Faces',f,...
        'FaceColor',varargin{ii},...
        'EdgeColor','none');
end
axis off;
background();
end

