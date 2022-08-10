function showcolor(varargin)
for ii = 1:nargin
    dx = ii-1;
    v = [dx,0;dx+1,0;dx+1,1;dx,1];
    f = [1,2,3,4];
    patch('Vertices',v,'Faces',f,...
        'FaceColor',varargin{ii},...
        'EdgeColor','none');
end
axis off;
background();
end

