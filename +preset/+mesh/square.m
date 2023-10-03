function msh = square(varargin)

    p = inputParser;
    addOptional(p,'n',50);
    addOptional(p,'w',1);
    parse(p,varargin{:});

    N = p.Results.n;
    W = p.Results.w;

    try
        sdf = sRectangle(W);
        msh = Mesh(sdf,'NElem',N);
        msh = msh.generate();
    end

end