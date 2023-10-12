function fem = beam(varargin)

p = inputParser;
addOptional(p,'n',50);
addOptional(p,'w',100);
addOptional(p,'h',5);
addOptional(p,'dt',1/300);
parse(p,varargin{:});

% try
    sdf = sRectangle(0, p.Results.w, 0,  p.Results.h);
    msh = Mesh(sdf,'Quads',[p.Results.n,3]);
    % msh = Mesh(sdf,'NElem',[p.Results.n]);
    msh = msh.generate();
    
    fem = Fem(msh,'TimeStep',p.Results.dt);
    fem = fem.addMaterial(NeoHookean(0.05, 0.45));
    fem = fem.addSupport('left', [1, 1]);
    fem = fem.addGravity();
    
    fem.solver.MaxIteration = 4;
% end

end