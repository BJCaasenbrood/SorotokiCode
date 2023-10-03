function msh = circle(varargin)

    p = inputParser;
    addOptional(p,'n',50);
    addOptional(p,'r',1);
    parse(p,varargin{:});

    N = p.Results.n;
    R = p.Results.r;

    try
        sdf = sCircle(R);
        msh = Mesh(sdf,'NElem',N);
        msh = msh.generate();
    end

end