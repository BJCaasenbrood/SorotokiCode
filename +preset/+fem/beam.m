function fem = beam(varargin)

p = inputParser;
addOptional(p,'n',50);
addOptional(p,'w',100);
addOptional(p,'h',5);
addOptional(p,'dt',1/300);
parse(p,varargin{:});

try
    sdf = sRectangle(0, p.Results.w, 0,  p.Results.h);
    msh = Mesh(sdf,'Quads',[p.Results.n,1]);
    % msh = Mesh(sdf,'NElem',[p.Results.n]);
    msh = msh.generate();
    
    fem = Fem(msh,'TimeStep',p.Results.dt);
    fem = fem.addMaterial(NeoHookean(0.1, 0.45));
    fem = fem.addSupport('left', [1, 1]);
    fem = fem.addGravity();

    % fem = fem.addContact(sCircle(10,[70,-40]));
    
    fem.solver.MaxIteration = 4;
    %fem.options.Display = [];
    %fem.solver.isLog  = false;

%    fem = solveDynamicFem(fem);
end

end