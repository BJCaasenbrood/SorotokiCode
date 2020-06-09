classdef Fem < handle

    properties (Access = public)
        Mesh;
        Material;
        Density;
        Node;
        Element;
        NNode;
        NElem;
        Dim;
        BdBox;
        Center;
        Topology;
        Log;
    end
    
    properties (Access = private)
        Node0;
        Support;
        Load;
        Spring;
        Output;
        Pressure;
        PressureCell;
        Contact;
        Contraction;
        ElemNDof;
        ElemSort;
        Normal;
        Edge;
        
        VonMisesNodal; 
        sxxNodal; syyNodal; sxyNodal;
        exxNodal; eyyNodal; exyNodal;
        fxNodal;  fyNodal;
        SPxyNodal;
        ElemMat;
        Rho = 1e12;
        Zeta = 0.1;
        fInternal = rand(1)*1e-8;
        fExternal = rand(1)*1e-8;
        Stiffness; TangentStiffness;
        Residual = Inf;
        
        Time = 0;
        TimeEnd = 1.0;
        TimeStep = 0.1;
        TimeStep0 = 0.1;
        TimeDelta = 0;
        TimeStepMin = 1e-3;
        EndIncrement;
        LoadingFactor;
        SigmoidFactor;
        
        ResidualNorm = 1e-3;
        StressNorm = 1e-9;
        DisplaceNorm = 1e-9;
        Objective = 1e-4;
        Constraint = 1e-4;
        
        IterationMMA = 1;
        Iteration = 1;
        Increment = 1;
        Divergence = 0;
        ShapeFnc;
        U = 0;
        Utmp = 0;
        
        MaxIteration = 50;
        MaxIterationMMA = 50;
        SolverStart = false;
        SolverStartMMA = false;
        SolverPlot = true;
        Assemble = false;
        Convergence = false; 
        Nonlinear = true;
        BisectLimit = 5;
        BisectCounter = 1;
        AssembledSystem = false;
        PrescribedDisplacement = false;
        VolumetricPressure = false;
        PressureLoad = 0;
        Type = 'PlaneStrain'
        Linestyle = '-';
        Colormap = turbo;
        I3 = eye(3); O3 = zeros(3);
        i; j; m; fi; t; e; c; s; p; v; l; k; fb; ed; fb0; ft;
        
        SolverResidual = 1e7;
        SolverVonMises = 1e7;
        SolverDisplace = 1e7;
        
        InformationBoolean = false;
        
        Movie = false;
        MovieStart = false;
        MovieAxis; MovieCAxis;
        
        VolumeInfill = 0.3;
        Ersatz = 1e-3; 
        Penal = 1;
        PenalMax = 4;
        PenalStep = 10;
        ChangeMax = 0.03;
        Beta = 1.5;
        FilterRadius = 0.5;
        SpatialFilter;
        OptFactor = 1;
        OptimizationSolver = 'MMA';
        MaterialInterpolation = 'SIMP';
        OptimizationProblem = 'Compliance';
        xold1; xold2; upp; low; fnorm;
        zMin = 0; zMax = 1;
        OutputVector; Change; dFdE;
        Periodic;
               
        VoidTolerance = 0.05*6; 
        
        Crop;
        Reflect;
        NumGridX; NumGridY;
        Grid; Grideval;
        Obj = rand(1)*1e-8;
        Con = rand(1)*1e-8;
        
        TopologyGridX; TopologyGridY; TopologyGridZ;
        TopologyGapFill;
        Repeat; ReflectionPlane; CellRepetion;
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Fem(Mesh,varargin) 
    
    obj.Mesh     = Mesh;
    obj.Node     = Mesh.get('Node');
    obj.Node0    = Mesh.get('Node');
    obj.NNode    = Mesh.get('NNode');
    obj.Dim      = Mesh.get('Dim');
    obj.Element  = Mesh.get('Element');
    obj.NElem    = Mesh.get('NElem');
    obj.Center   = Mesh.get('Center');
    obj.BdBox    = Mesh.get('BdBox');
    obj.Density  = ones(obj.NElem,1);
    obj.Residual = zeros(obj.Dim*obj.NNode,1);
    obj.Utmp     = zeros(obj.Dim*obj.NNode,1);
    obj.SigmoidFactor  = 0;
    obj.Load     = zeros([1,0,0]);
       
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    obj.TimeStep0 = obj.TimeStep;
    obj = SetupFiniteElement(obj);
end

%---------------------------------------------------------------------- get     
function varargout = get(Fem,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Fem.(varargin{ii});
        end
    else
        varargout = Fem.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function Fem = set(Fem,varargin)
    for ii = 1:2:length(varargin)
        Fem.(varargin{ii}) = varargin{ii+1};
        if strcmp(varargin{ii},'TimeStep')
            Fem.TimeStep0 = varargin{ii+1}; 
        end
        if strcmp(varargin{ii},'FilterRadius')
            Fem.SpatialFilter = GenerateRadialFilter(Fem,varargin{ii+1});
        end
    end
    
end

%--------------------------------------------------------------------- show
function h = show(Fem,varargin)
    
if nargin<2, Request = '0'; 
else, Request = varargin{1}; end

S = 'interp'; 
V = Fem.Node;
flag = 0;
colormap(Fem.Colormap);

switch(Request)
    case('0'),   Z = zeros(Fem.NNode,1);
    case('SDF'), Z = Fem.Mesh.SDF(Fem.Mesh.Node); Z = Z(:,end);
    case('Svm'), Z = Fem.VonMisesNodal;
    case('Sxx'), Z = Fem.sxxNodal;
    case('Syy'), Z = Fem.syyNodal;
    case('Sxy'), Z = Fem.sxyNodal;
    case('Exx'), Z = Fem.exxNodal;
    case('Eyy'), Z = Fem.eyyNodal;
    case('Exy'), Z = Fem.exyNodal;
    case('Fx'),  Z = Fem.fxNodal;
    case('Fy'),  Z = Fem.fyNodal;
    case('Fin'), [~,~,Z] = DisplacementField(Fem,Fem.fInternal);
    case('Fex'), [~,~,Z] = DisplacementField(Fem,Fem.fExternal);
    case('Un'),  [~,~,Z] = DisplacementField(Fem,Fem.Utmp);
    case('Ux'),  [Z,~,~] = DisplacementField(Fem,Fem.Utmp);
    case('Uy'),  [~,Z,~] = DisplacementField(Fem,Fem.Utmp);
    case('E'),   [~,~,Z] = MaterialField(Fem); S = 'flat'; 
    V = Fem.Node0; colormap(noir(-1)); background('w');
    case('E+'),  [~,~,Z] = MaterialField(Fem); S = 'flat'; 
    colormap(noir(-1)); background('w');
    otherwise;   flag = 1; Z = 0;
end

if length(Z) == Fem.NNode, Z = Smoothing(Fem,Z,1); end

if length(varargin) == 2 && flag == 0
    handle = varargin{2}; 
    set(handle{2},'FaceVertexCData',Z); 
    set(handle{2},'Vertices',V); 
    return; 
end

if flag == 0
cla; axis equal;     
axis off; hold on; h{3} = [];

FaceMatrix  = Fem.Mesh.get('ElemMat');
BoundMatrix = Fem.Mesh.get('Boundary');

if Fem.Dim > 2 && (strcmp(Request,'E') || strcmp(Request,'E+'))
    Z(Z<Fem.VoidTolerance,:) = nan;
    S = 'flat';
    %Alp = Fem.Mesh.get('ElementToFace')*Z;
    Z = Fem.Mesh.get('ElementToFace')*Z;
else
    Alp = ones(Fem.NElem,1);
end


h{1} = patch('Faces',BoundMatrix,'Vertices',Fem.Node0,...
    'LineStyle','-','Linewidth',1,'EdgeColor','k');

h{2} = patch('Faces',FaceMatrix,'Vertices',V,...
    'FaceVertexCData',Z,'Facecolor',S,'LineStyle',Fem.Linestyle,...
    'Linewidth',1.5,'FaceAlpha',1,'EdgeColor','k');

h{3} = patch('Faces',BoundMatrix,'Vertices',V,...
    'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

if Fem.Dim == 3, view(30,10); end

if Fem.VolumetricPressure
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density); 
    if strcmp(Request,'E+'), Pc = computeCentroid(Fem);
    else, Pc = Fem.Mesh.get('Center'); 
    end
    h{5} = plot(Pc(id,1),Pc(id,2),'.','Color','k');
end

if ~isempty(Fem.Contact)
    SDF = Fem.Contact{1};
    Move = Fem.Contact{2};
    BD = 2*Fem.BdBox;
    [px,py] = meshgrid(linspace(BD(1),BD(2),50),linspace(BD(3),BD(4),50));
    Y = [px(:),py(:)];
    beta = 0.975*Fem.Time;
    Y(:,1) = Y(:,1) - beta*Move(1);
    Y(:,2) = Y(:,2) - beta*Move(2);
    d = SDF(Y);
    D = reshape(d(:,end),[50,50]);
    h{5} = contourf(px,py,-D,[1e-3 1e-3],'linewidth',1.5,...
        'FaceColor',0.85*gitpage);
end
end

if flag == 1
switch(Request)
case('Residual'),   semilogy(Fem.SolverResidual,'linewidth',1);
case('Stress'),     semilogy(Fem.SolverVonMises,'linewidth',1);
case('Displace'),   semilogy(Fem.SolverDisplace,'linewidth',1);
case('Objective'),  plot(Fem.Objective,'linewidth',1);
case('Constraint'), plot(Fem.Constraint,'linewidth',1);
otherwise,          flag = 2;
end
end

if flag == 2
    switch(Request)
    case('BC')
    if ~isempty(Fem.Support),  Supp = Fem.Support; else, Supp = []; end
    if ~isempty(Fem.Load),     Forc = Fem.Load; else, Forc = []; end
    if ~isempty(Fem.Spring),   Sprg = Fem.Spring(:,1); else, Sprg = []; end
    if ~isempty(Fem.Output),   Out = Fem.Output; else, Out = []; end
    if ~isempty(Fem.Pressure), Pres = Fem.Pressure(:,1); else, Pres = []; end
    
    hold on;
    for ii = 1:size(Supp,1)
        id = Supp(ii,1);
        symbol = '\Delta';
        if Fem.Dim == 2
        plot(V(id,1),V(id,2),'marker','d','markersize',5,'Color',col(1));
        else
        plot3(V(id,1),V(id,2),V(id,3),'marker','d',...
            'markersize',5,'Color',col(1));
        end
    end
    
    for ii = 1:size(Forc,1)
        id = Forc(ii,1);
        if ispos(Forc(ii,2)), symbol = '\bf \rightarrow'; end
        if isneg(Forc(ii,2)), symbol = '\bf \leftarrow'; end
        if ispos(Forc(ii,3)), symbol = '\bf \uparrow'; end
        if isneg(Forc(ii,3)), symbol = '\bf \downarrow'; end
        vecmag = 0.15*sum(abs(axis))/4;
    end
    
    plotvector(V(Forc(:,1),:),Forc(:,2:end)*vecmag);
    
    for ii = 1:size(Out,1)
        id = Out(ii);
        if ispos(Out(ii,2)), symbol = '\bf \rightarrow'; end
        if isneg(Out(ii,2)), symbol = '\bf \leftarrow'; end
        if ispos(Out(ii,3)), symbol = '\bf \uparrow'; end
        if isneg(Out(ii,3)), symbol = '\bf \downarrow'; end
        text(V(id,1)-.1,V(id,2)-.1,symbol,'fontsize',10,'Color','g');
    end
    
    for ii = 1:length(Sprg)
        id = Sprg(ii);
        text(V(id,1),V(id,2),symbol,'fontname',font,'fontsize',10,'Color','g');
    end
    
        
    for ii = 1:length(Pres)
        id = Pres(ii);
        text(V(id,1),V(id,2)-.1,symbol,'fontsize',10,'Color','g');
    end
    
    end
    
    flag = 3;
end

if flag == 3
    switch(Request)
    case('ISO')
        clf; former(Fem,10);
        if nargin < 3, ISOVALUE = 0.2; else, ISOVALUE = varargin{2};
        end
        showISO(Fem,ISOVALUE,0.5);
        %showISOPlus(Fem,ISOVALUE,0.5);
    case('ISO+')
        clf; former(Fem,10);
        if nargin < 3, ISOVALUE = 0.2; else, ISOVALUE = varargin{2};
        end
        %showISO(Fem,ISOVALUE,0.5);
        showISOPlus(Fem,ISOVALUE,0.5);
    end    
end

if Fem.Movie %&& strcmp(Request,'ISO')
    background(metropolis);
    if Fem.MovieStart == false
       Fem.MovieStart = true;
       if ~Fem.SolverStartMMA, Name = 'fem'; else, Name = 'topo'; end
       MovieMaker(Fem,Name,'Start');
    else
       MovieMaker(Fem);
    end
end

end

%-------------------------------------------------------------------- reset
function Fem = reset(Fem,Request)
   switch(Request) 
       case('fem'), Fem.Iteration = 1; Fem.Increment = 1; 
           Fem.Node = Fem.Node0; Fem.Utmp = [];
           Fem.Spring = [];
       case('opt'), Fem.Iteration = 1; Fem.Increment = 1; 
           Fem.IterationMMA = 1; Fem.Node = Fem.Node0; 
           Fem.Utmp = []; Fem.Change = []; Fem.Change;
           Fem.Objective = []; Fem.Constraint = [];
   end
    
end

%-------------------------------------------------------------------- solve
function Fem = solve(Fem)
    
Fem.Log           = [];    
Fem.TimeDelta     = Fem.TimeEnd;
Fem.BisectCounter = 1;
Fem.Increment     = 1;
Fem.Iteration     = 1;
Fem.Time          = 0;
Fem.TimeDelta     = 0;
Fem.EndIncrement  = false;
Fem.Convergence   = true;
Fem.TimeStep      = Fem.TimeStep0 + 1e-6;
Fem.Density       = clamp(Fem.Density,Fem.Ersatz,1);
Fem.Utmp          = zeros(Fem.Dim*Fem.NNode,1);
Fem.Node          = Fem.Node0;

showInformation(Fem,'NonlinearFem');

if Fem.SolverPlot || ~Fem.SolverStartMMA
    figure(101); 
    Fem.show('0');
    background(metropolis);
end

while true 
    
    [Fem,Terminate] = SolverTimer(Fem);
    if Terminate, break; end
        
    flag = 0;
    Singular = false;
    Fem.Divergence = 0;
        
    % load increment
    Fem.LoadingFactor = sigmoid(Fem.Time/Fem.TimeEnd,Fem.SigmoidFactor);
    
    % reset iteration
    Fem.Iteration = 1;
        
    while flag == 0
        
        % assemble global system matrices
        Fem = AssembleGlobalSystem(Fem);
        
        % assemble force/boundary conditions
        Fem = AssembleBoundaryConditions(Fem);
        
        % get free DOFs
        FreeDofs = GetFreeDofs(Fem);
        
        % compute residual and tangent stiffness
        if Fem.Nonlinear
            A = Fem.TangentStiffness(FreeDofs,FreeDofs);
            B = Fem.Residual(FreeDofs);
        elseif ~Fem.Nonlinear
            A = Fem.Stiffness(FreeDofs,FreeDofs);
            B = sparse(Fem.fExternal(FreeDofs));
        end
        
        Delta = Fem.Utmp;

        if rcond(full(A)) >= 1e-10, DeltaU = A\B;
        else, Singular = true; DeltaU = Fem.Utmp(FreeDofs)*0;
        end
            
        if Fem.Nonlinear, Delta(FreeDofs,1)=Delta(FreeDofs,1)-DeltaU(:,1);
        else, Delta(FreeDofs,1) = DeltaU(:,1); B = Fem.ResidualNorm; end
        
        % evaluate convergence criteria
        Fem.SolverResidual = vappend(Fem.SolverResidual,...
            max(abs(Fem.Residual(FreeDofs))));
        
        Fem.SolverVonMises = vappend(Fem.SolverVonMises,...
            max(Fem.s(:,1)));
        
        Fem.SolverDisplace = vappend(Fem.SolverDisplace,...
            max(DeltaU(:,1)));
        
        % check convergence
        [flag,Fem] = CheckConvergence(Fem,Singular); 
        
        % add iteration
        Fem.Iteration = Fem.Iteration + 1;
        
        % update mesh 
        if flag == 0 && Fem.Nonlinear
            Fem.Utmp = Delta;
        elseif flag == 1 && ~Fem.Nonlinear
            Fem.Utmp = Delta;
            Fem.Node = UpdateNode(Fem,Delta);
        elseif flag == 1 && Fem.Nonlinear
            Fem.Utmp = Delta;
            Fem.Node = UpdateNode(Fem,Delta);
        end
    end 
    
    if ~Fem.SolverStartMMA
        Fem = AssembleGlobalSystem(Fem, true);
        
        Fem.VonMisesNodal = full(sparse(Fem.l,1,Fem.s(:,1))./sparse(Fem.l,1,Fem.v));
        Fem.sxxNodal = full(sparse(Fem.l,1,Fem.s(:,2))./sparse(Fem.l,1,Fem.v));
        Fem.syyNodal = full(sparse(Fem.l,1,Fem.s(:,3))./sparse(Fem.l,1,Fem.v));
        Fem.sxyNodal = full(sparse(Fem.l,1,Fem.s(:,4))./sparse(Fem.l,1,Fem.v));
        Fem.exxNodal = full(sparse(Fem.l,1,Fem.p(:,1))./sparse(Fem.l,1,Fem.v));
        Fem.eyyNodal = full(sparse(Fem.l,1,Fem.p(:,2))./sparse(Fem.l,1,Fem.v));
        Fem.exyNodal = full(sparse(Fem.l,1,Fem.p(:,3))./sparse(Fem.l,1,Fem.v));
        force = full(sparse(Fem.i,1,Fem.fi));
        Fem.fxNodal = force(2*(1:Fem.NNode)-1);
        Fem.fyNodal = force(2*(1:Fem.NNode));
    end
    
    if ~isempty(Fem.Output) && ~Fem.SolverStartMMA
       [ux,uy,un] = DisplacementField(Fem,Fem.Utmp);
       [fx,fy,fn] = DisplacementField(Fem,Fem.fInternal);
       idNodes = Fem.Output(:,1);
       if isempty(Fem.Log)
           Fem.Log = cell(10,2);
           Fem.Log{1,1} = 'ux'; Fem.Log{1,2} = ux(idNodes); 
           Fem.Log{2,1} = 'uy'; Fem.Log{2,2} = uy(idNodes);
           Fem.Log{3,1} = 'un'; Fem.Log{3,2} = un(idNodes);
           Fem.Log{4,1} = 'fx'; Fem.Log{4,2} = fx(idNodes); 
           Fem.Log{5,1} = 'fy'; Fem.Log{5,2} = fy(idNodes);
           Fem.Log{6,1} = 'fy'; Fem.Log{6,2} = fn(idNodes);
           Fem.Log{7,1} = 'Svm'; Fem.Log{7,2} = Fem.VonMisesNodal(idNodes);
           Fem.Log{8,1} = 'Sxx'; Fem.Log{8,2} = Fem.sxxNodal(idNodes);
           Fem.Log{9,1} = 'Syy'; Fem.Log{9,2} = Fem.syyNodal(idNodes);
           Fem.Log{10,1} = 'Sxy'; Fem.Log{10,2} = Fem.sxyNodal(idNodes);
       else
           Fem.Log{1,2} = vappend(Fem.Log{1,2},ux(idNodes));
           Fem.Log{2,2} = vappend(Fem.Log{2,2},uy(idNodes));
           Fem.Log{3,2} = vappend(Fem.Log{3,2},un(idNodes));
           Fem.Log{4,2} = vappend(Fem.Log{4,2},fx(idNodes));
           Fem.Log{5,2} = vappend(Fem.Log{5,2},fy(idNodes));
           Fem.Log{6,2} = vappend(Fem.Log{6,2},fn(idNodes));
           Fem.Log{7,2} = vappend(Fem.Log{7,2},Fem.VonMisesNodal(idNodes));
           Fem.Log{8,2} = vappend(Fem.Log{8,2},Fem.sxxNodal(idNodes));
           Fem.Log{9,2} = vappend(Fem.Log{9,2},Fem.syyNodal(idNodes));
           Fem.Log{10,2} = vappend(Fem.Log{10,2},Fem.sxyNodal(idNodes));
       end
    end
    
    if Fem.SolverPlot || ~Fem.SolverStartMMA
        Fem.show('Svm'); 
        drawnow;
    end
    
    if ~Fem.Nonlinear, break; end
end

end

%------------------------------------------------------ optimization solver
function Fem = optimize(Fem)
 
showInformation(Fem,'TopologyOptimization');
    
Fem.SpatialFilter  = GenerateRadialFilter(Fem,Fem.FilterRadius);
Fem.SolverPlot     = false;
Fem.IterationMMA   = 0;
Fem.SolverStartMMA = true;
flag               = true;
Visual             = 'ISO';

% draw initial visual
Fem.show(Visual); drawnow;

while flag

    % update increment
    Fem.IterationMMA = Fem.IterationMMA + 1;
    
    % update material field
    [~,dEdy,~,dVdy] = MaterialField(Fem);
     
    % set load factor
    Fem.OptFactor = Fem.IterationMMA/Fem.MaxIterationMMA;
   
    % solve nonlinear finite elements
    Fem = Fem.solve();
    
    % evaluate objective and constraints
    [f,dfdE,dfdV] = ObjectiveFunction(Fem);
    [g,dgdE,dgdV] = ConstraintFunction(Fem);
    
    % compute design sensitivities
    dfdz = Fem.SpatialFilter'*(dEdy.*dfdE + dVdy.*dfdV);
    dgdz = Fem.SpatialFilter'*(dEdy.*dgdE + dVdy.*dgdV);
    
    %compute design variable
    if strcmp(Fem.OptimizationSolver,'MMA')
       [Fem,ZNew] = UpdateSchemeMMA(Fem,f,dfdz,g,dgdz);
    elseif strcmp(Fem.OptimizationSolver,'OC')
        [Fem,ZNew] = UpdateSchemeOC(Fem,dfdz,g,dgdz);
    end
    
    % determine material change
    Fem.Change  = clamp(ZNew - Fem.Density,-Fem.ChangeMax,Fem.ChangeMax);
    Fem.Density = Fem.Density + Fem.Change;
    
    % evaluate fitness
    Fem.Objective  = vappend(Fem.Objective,f);       
    Fem.Constraint = vappend(Fem.Constraint,g);
    
    % check convergence
    [flag,Fem] = CheckConvergenceOpt(Fem);
    
    % draw visual
    Fem.show(Visual); drawnow;
end

Fem.SolverStartMMA = false;

end

%----------------------------------------------------- form smooth topology
function Fem = former(Fem, Thickness)
Res = 300;    
Layers = 40;
Patch = 5;

if nargin < 2, Thickness = 0.2*(Fem.BdBox(2)-Fem.BdBox(1)); end
verts = Fem.Node0; 

V = Fem.Mesh.get('NodeToFace')*Fem.SpatialFilter*Fem.Density;

if Fem.VolumetricPressure
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
    Pc = Fem.Mesh.Center;
    S = zeros(Fem.NElem,1); 
    d = DistancePointSet(Fem,Pc,Pc(id,:),0.5*Fem.FilterRadius);
    ide = unique(d(:,2));
    rhotmp =  Fem.SpatialFilter*Fem.Density;
    S(ide,1) = (rhotmp(ide)>Fem.VoidTolerance)*0.5;
    Fill = zeros(Fem.NElem,1); Fill(id) = 0.5;
    E = Fem.Mesh.get('NodeToFace')*...
   Fem.SpatialFilter*(S + Fem.Density*0.01);
    V2 = clamp(Fem.Mesh.get('NodeToFace')*...
   Fem.SpatialFilter*(Fill + Fem.Density-0.1),0,1);
end

if isempty(Fem.Crop), B = Fem.BdBox; 
else, B = Fem.Crop;
end

x = verts(:,1); y = verts(:,2);
dx = 0.01*(B(2) - B(1));
dy = 0.01*(B(4) - B(3));
xq = linspace(B(1)+dx,B(2)-dx,Res); 
yq = linspace(B(3)+dy,B(4)-dy,Res);

[xxq,yyq] = meshgrid(xq,yq);
P = griddata(x,y,V,xxq,yyq);
vrt = [xxq(:),yyq(:)];
Dist = Fem.Mesh.SDF(vrt); Dist = Dist(:,end);
tol = 0.1*sqrt((B(2)-B(1))*(B(4)-B(3))/length(vrt));
P = P(:); P(Dist > tol) = 0; 
P = (reshape(P,Res,Res)); 
V = cat(3,xxq,yyq,P);

if Fem.VolumetricPressure
    P = griddata(x,y,V2,xxq,yyq);
    Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); Dist = Dist(:,end);
    P = P(:); P(Dist > tol) = 0; P(isnan(P))=0;
    P = (reshape(P,Res,Res));
    
    GapFill = cat(3,xxq,yyq,P);
    
    P = griddata(x,y,E,xxq,yyq);
    Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); Dist = Dist(:,end);
    P = P(:);  P(isnan(P))=0;
    P = (reshape(P,Res,Res));
    
    EdgeFill = cat(3,xxq,yyq,P);
end

SDF0 = V(:,:,3);
SDF0 = GaussianFilter(SDF0,5);
SDF = repmat(SDF0,[1 1 Layers]);

if Fem.VolumetricPressure
    V = GapFill;
    SDF2 = V(:,:,3);
    SDF2 = GaussianFilter(SDF2,10);
    ZFiller = repmat(SDF2,[1 1 Patch]);
    for ii = 1:Patch
        ZFiller(:,:,ii) = lerp(SDF2+0.01*ii,SDF0,cos((1/Patch)*pi));
    end
    %SDF(:,:,2:Patch+1) = ZFiller;
end

% zero padding
SDF(:,:,1) = 0; SDF(:,:,end) = 0;
SDF(1,:,:) = 0; SDF(end,:,:) = 0;
SDF(:,1,:) = 0; SDF(:,end,:) = 0;

Fem.Topology = SDF;

zq = linspace(0,Thickness,Layers);
[X,Y,Z] = meshgrid(xq,yq,zq);

Fem.TopologyGridX = X; 
Fem.TopologyGridY = Y; 
Fem.TopologyGridZ = Z; 

end

%--------------------------------------------------------- show iso-surface
function Fem = showISO(Fem,varargin)
V = Fem.Topology;
X = Fem.TopologyGridX; 
Y = Fem.TopologyGridY; 
a = size(V);

if nargin < 3, depth = floor(a(3)/2);
else, depth = a(3)*varargin{2}; end

X = X(:,:,depth);
Y = Y(:,:,depth);
SDF0 = V(:,:,depth);

SDF = GenerateCell(Fem,SDF0);
if ~isempty(Fem.CellRepetion)
    Rpt = Fem.CellRepetion;
    SDF = repmat(SDF,Rpt(2),Rpt(1));
else, Rpt = [1,1];
end

if ~isempty(Fem.ReflectionPlane)
Rp = Fem.ReflectionPlane;
if Rp(1) == 1
    Xf = flip(X) + (max(max(max(X))) - min(min(min(X))));
    X = vertcat(Xf,X);
end
if Rp(1) == -1
    Xf = flip(X) - (max(max(max(X))) - min(min(min(X))));
    X = horzcat(X,Xf);
end
if Rp(2) == -1
    Yf = flip(Y) - (max(max(max(Y))) - min(min(min(Y))));
    Y = horzcat(Y,Yf);
end
if Rp(2) == 1
    Yf = flip(Y) + (max(max(max(Y))) - min(min(min(Y))));
    Y = horzcat(Yf,Y);
end
end

scaleX = 1; scaleY = 1;
if ~isempty(Fem.Repeat)
instr = Fem.Repeat;
SDF0 = SDF;
for ii = 1:length(instr)
    if instr(ii) == 1
      scaleX = scaleX + 1;
      SDF = cat(2,SDF,SDF0);
    end
    if instr(ii) == 2
      scaleY = scaleY*2;
      SDF = cat(1,SDF,SDF);
    end
end
end

cla;
Uxx = scaleX*Rpt(1)*[min(min(min(X))) max(max(max(X)))];
Uyy = scaleY*Rpt(2)*[min(min(min(Y))) max(max(max(Y)))];
I = GaussianFilter((SDF >= varargin{1})*255,3);
[XX,YY] = meshgrid(linspace(Uxx(1),Uxx(2),size(I,1)),...
                   linspace(Uyy(1),Uyy(2),size(I,2)));

image(rescale(Uxx),...
    ((max(Uyy)-min(Uyy))/(max(Uxx) - min(Uxx)))*rescale(Uyy),I);

axis equal; axis off; 
colormap(noir(-1)); 
caxis([0 1]);
background('w');
end

%------------------------------------------------ show deformed iso-surface
function Fem = showISOPlus(Fem,varargin)
V = Fem.Topology;
X = Fem.TopologyGridX; 
Y = Fem.TopologyGridY; 
a = size(V);

if nargin < 3, depth = floor(a(3)/2);
else, depth = a(3)*varargin{2}; end

X = X(:,:,depth);
Y = Y(:,:,depth);
SDF0 = V(:,:,depth);

SDF = GenerateCell(Fem,SDF0);
if ~isempty(Fem.CellRepetion)
    Rpt = Fem.CellRepetion;
    SDF = repmat(SDF,Rpt(2),Rpt(1));
else, Rpt = [1,1];
end

if ~isempty(Fem.ReflectionPlane)
Rp = Fem.ReflectionPlane;
if Rp(1) == 1
    Xf = flip(X) + (max(max(max(X))) - min(min(min(X))));
    X = vertcat(Xf,X);
end
if Rp(1) == -1
    Xf = flip(X) - (max(max(max(X))) - min(min(min(X))));
    X = horzcat(X,Xf);
end
if Rp(2) == -1
    Yf = flip(Y) - (max(max(max(Y))) - min(min(min(Y))));
    Y = horzcat(Y,Yf);
end
if Rp(2) == 1
    Yf = flip(Y) + (max(max(max(Y))) - min(min(min(Y))));
    Y = horzcat(Yf,Y);
end
end

scaleX = 1; scaleY = 1;
if ~isempty(Fem.Repeat)
instr = Fem.Repeat;
SDF0 = SDF;
for ii = 1:length(instr)
    if instr(ii) == 1
      scaleX = scaleX + 1;
      SDF = cat(2,SDF,SDF0);
    end
end
end

cla;
Uxx = scaleX*Rpt(1)*[min(min(min(X))) max(max(max(X)))];
Uyy = scaleY*Rpt(2)*[min(min(min(Y))) max(max(max(Y)))];
I = GaussianFilter((SDF >= varargin{1})*255,3);

% Uyy = ((max(Uyy)-min(Uyy))/(max(Uxx) - min(Uxx)))*rescale(Uyy);             
% Uxx = rescale(Uxx); 
% 
% %image(Uxx,Uyy,I);
% 
% [XX,YY] = meshgrid(linspace(Uxx(1),Uxx(2),size(I,1)),...
%                    linspace(Uyy(1),Uyy(2),size(I,2)));

cla; axis equal;     
axis off; hold on;
%if length(Z) ~= Fem.NNode, T=Fem.Mesh.get('NodeToFace'); Z=T*Z; end
Pc = computeCentroid(Fem);
[VV0, f0] = QuickTriangulation(Fem,Fem.Mesh.get('Center'),Fem.Node0,Fem.Element);
FF0 = GenerateElementalMatrix(f0,length(f0));
[VV, f] = QuickTriangulation(Fem,Pc,Fem.Node,Fem.Element);
FF = GenerateElementalMatrix(f,length(f));

VV = [VV, VV(:,1)*0];
%VT = [XX(:), YY(:)];

Uyy = rescale(VV0(:,2));             
Uxx = rescale(VV0(:,1)); 

patcht(FF,VV,FF0,[Uyy, Uxx],rot90(I'));

patch('Faces',Fem.Mesh.get('Boundary'),'Vertices',Fem.Node,...
    'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

axis equal; axis off; 
colormap(metro(-1)); 
caxis([0 1]);
background('w');

end

%-------------------------------------------------- show STL reconstruction
function obj = showSTL(Fem,value)
    
V = Fem.Topology;
X = Fem.TopologyGridX; 
Y = Fem.TopologyGridY; 
Z = Fem.TopologyGridZ; 

[faces,vertices] = MarchingCubes(X,Y,Z,V,value);

figure(101);
cla;
obj = Gmodel(vertices,faces);
obj.Texture = grey;

obj = obj.bake();
obj = Blender(obj,'Rotate',{'x',90});
obj.render();
%
%         contourf(Z,[value value],...
%             'linestyle','-','EdgeColor','k','Linewidth',2), 
%         
%         patch('Faces',Fem.Mesh.get('Boundary'),'Vertices',vert,...
%          'LineStyle','-','Linewidth',2,'EdgeColor','k');
%      
axis equal
axis off;
%colormap(inferno);
%caxis([-1 1.5]);
%     end

end

%----------------------------------------------------------- initial design
function Fem = initialTopology(Fem,varargin)

    if strcmp(varargin{1},'Equidistance')
        Fem.Density = InitialDesign(Fem,{varargin{2},varargin{3}});
    elseif strcmp(varargin{1},'Hole')
        Pc = Fem.Mesh.get('Center');
        d = DistancePointSet(Fem,varargin{2},Pc,varargin{3});
        Z = ones(Fem.NElem,1); Z(d(:,1)) = 0;
        Fem.Density = Z;
    else
        Fem.Density = InitialDesign(Fem,{[1,1],1});
    end
    
    Fem.Density = Fem.SpatialFilter*(1.25*Fem.Density*Fem.VolumeInfill);
end

%--------------------------------------------------------------- find nodes
function NodeList = FindNodes(Fem,varargin)
    if  nargin > 2
        NodeList = FindNode(Fem.Node,varargin{1:end});
    else
        NodeList = FindNode(Fem.Node,varargin{1:end});
    end
end

%------------------------------------------------------------ find elements
function ElementList = FindElements(Fem,varargin)
    ElementList = FindNode(Fem.Mesh.Center,varargin{1:end});
end

%---------------------------------------------------------- add constraints
function Fem = AddConstraint(Fem,varargin)
    
for ii = 1:3:length(varargin)
  if size(varargin{ii+2},2) == 2 || size(varargin{ii+2},2) == 3
      if strcmp(varargin{ii},'PressureCell')
         Fem.VolumetricPressure = true; 
         Fem.(varargin{ii}) = [varargin{ii+1},repmat(varargin{ii+2},...
          [length(varargin{ii+1}),1])];
      elseif strcmp(varargin{ii},'Contact')
          Fem.(varargin{ii}) = {varargin{ii+1},varargin{ii+2}};
      else 
          BC = [varargin{ii+1},repmat(transpose(varargin{ii+2}(:)),...
          [length(varargin{ii+1}),1])];
            Fem.(varargin{ii}) = [Fem.(varargin{ii});BC];
      end
  else
      warning([varargin{ii}, ' has incorrect input'] );
  end
end
end

end
methods (Access = private)
%%///////////////////////////////////////////////////////// FINITE ELEMENTS
%------------------------------------------------ convert mesh to fem class
function Fem = SetupFiniteElement(Fem)

Fem.SpatialFilter         = GenerateRadialFilter(Fem,Fem.FilterRadius);
[~,~,Fem.Normal,Fem.Edge] = computeCentroid(Fem);

Fem.ElemNDof = Fem.Dim*cellfun(@length,Fem.Element);

if (~Fem.AssembledSystem && ~Fem.Nonlinear) || Fem.Nonlinear
Fem.i = zeros(sum(Fem.ElemNDof.^2),1);
Fem.j  = Fem.i; Fem.e  = Fem.i; Fem.fi = Fem.i; Fem.k  = Fem.i; 
Fem.m  = Fem.i; Fem.c  = Fem.i; Fem.t  = Fem.i; Fem.fb = Fem.i; 
Fem.ft = Fem.i; Fem.v  = Fem.l;
Fem.s  = zeros(Fem.NNode,6); 
Fem.s  = zeros(Fem.NNode,3);
Fem.l  = zeros(Fem.NNode,1); 
end

end

%------------------------------------ assemble global finite-element system 
function Fem = AssembleGlobalSystem(Fem,ForceBuild)
    
if nargin < 2, ForceBuild = false; end

% evaluate shape-functions at nodal locations
if (Fem.Iteration == 1 && Fem.Increment == 1)
    tmp = struct; tmp.Element = Fem.Element;
    switch(Fem.Mesh.Type)
        case('C2PX'), tmp = TabulateShapeFunctions(tmp);
        case('C2T3'), tmp = TabulateShapeFunctions(tmp);
        case('C2Q4'), tmp = TabulateShapeFunctions(tmp);
        case('C3H8'), tmp = TabulateShapeFunctionsC3H8(tmp);
        otherwise, tmp = TabulateShapeFunctions(tmp);
    end
    Fem.ShapeFnc = tmp.ShapeFnc;
end

if Fem.VolumetricPressure
    beta = Fem.OptFactor*Fem.LoadingFactor;
    dV = beta*Fem.PressureCell(:,2);
else, dV = 0; 
end

[E,~,V] = MaterialField(Fem);

index    = 0; 
subindex = 0;

if (~Fem.AssembledSystem && ~Fem.Nonlinear) || Fem.Nonlinear || ForceBuild
    
for el = 1:Fem.NElem
   
    NDof = Fem.ElemNDof(el);
    if Fem.Dim == 2
    eDof = reshape([2*Fem.Element{el}-1; 2*Fem.Element{el}],NDof,1);
    else
    eDof = reshape([3*Fem.Element{el}-2; 3*Fem.Element{el}-1;
                    3*Fem.Element{el}],NDof,1);
    end
    
    [Fe,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te] = ...
    Locals(Fem,Fem.Element{el},eDof,dV,V(el));
    
    ind1 = index+1:index+NDof^2;
    ind2 = index+1:index+NDof;
    ind3 = subindex+1:subindex+NDof/Fem.Dim;

    I = repmat(eDof,1,NDof); J = I';
    Fem.e(ind1)  = el;
    Fem.i(ind1)  = I(:);
    Fem.j(ind1)  = J(:);
    Fem.m(ind1)  = Me(:);
    Fem.c(ind1)  = Ce(:);
    Fem.k(ind1)  = Ke(:);
    Fem.t(ind1)  = Kte(:);
    Fem.fi(ind2) = Fe(:);
    Fem.ft(ind2) = Te(:);    
    
    Fem.s(ind3,1) = E(el)*Svme(:);
    Fem.s(ind3,2) = E(el)*SS(:,1);
    Fem.s(ind3,3) = E(el)*SS(:,2);
    Fem.s(ind3,4) = E(el)*SS(:,4);
    Fem.p(ind3,1) = E(el)*EE(:,1);
    Fem.p(ind3,2) = E(el)*EE(:,2);
    Fem.p(ind3,3) = E(el)*EE(:,4);
    Fem.v(ind3) = Qe(:);
    Fem.l(ind3) = Fem.Element{el}(:);
    
    index    = index + NDof^2;
    subindex = subindex + NDof/Fem.Dim;
end
end

Fem.AssembledSystem = true;

end

%------------------------ assemble prescribed boundary forces/displacements
function Fem = AssembleBoundaryConditions(Fem)
    
[E,~,~,~] = MaterialField(Fem);

% build spring/elastic forces
NSpring = size(Fem.Spring,1);
sp      = sparse(Fem.Dim*Fem.NNode,Fem.Dim);
if NSpring ~=0
    sp(2*Fem.Spring(1:NSpring)-1) = Fem.Spring(1:NSpring,2);
    sp(2*Fem.Spring(1:NSpring))   = Fem.Spring(1:NSpring,3);
end

spMat = spdiags(sp(:),0,Fem.Dim*Fem.NNode,Fem.Dim*Fem.NNode);

% build output vector of nodal displacements
NOutput = size(Fem.Output,1);
L       = sparse(Fem.Dim*Fem.NNode,1);
if NOutput ~=0
    L(2*Fem.Output(1:NOutput,1)-1) = Fem.Output(1:NOutput,2);
    L(2*Fem.Output(1:NOutput,1)) = Fem.Output(1:NOutput,3);
end
    
K   = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
Ktr = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
Fe  = sparse(Fem.i,1,E(Fem.e).*Fem.fi);
F   = sparse(Fem.Dim*Fem.NNode,1);

if Fem.Nonlinear, beta = Fem.OptFactor*Fem.LoadingFactor;
else, beta = 1;
end

if ~isempty(Fem.Load) && ~Fem.PrescribedDisplacement
NLoad = size(Fem.Load,1);
% if ~isempty(Fem.Mesh.get('NodeToFace')) 
% A = sign(Fem.Mesh.get('NodeToFace'));
% Gamma = A(Fem.Load(1:NLoad,1),:);
% W = sum(Gamma,2)/sum(sum(Gamma,2));
% else, W = 1.0;
W = 1.0;
%end

for ii = 1:Fem.Dim
F(Fem.Dim*Fem.Load(1:NLoad,1)+(ii-Fem.Dim),1) = beta*W.*Fem.Load(1:NLoad,1+ii);
end

end

if ~isempty(Fem.Contact)
    SDF = Fem.Contact{1};
    Move = Fem.Contact{2};
    Y = Fem.Node;
    Y(:,1) = Y(:,1) - beta*Move(1);
    Y(:,2) = Y(:,2) - beta*Move(2);
    d = SDF(Y); I = find((d(:,end))<0);   
    eps = 1e-5;
    n1 = (SDF(Y(I,:)+repmat([eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
    n2 = (SDF(Y(I,:)+repmat([0,eps],size(Y(I,:),1),1))-d(I,end))/eps;
    Ux = -1.95*d(I,end).*n1(:,end);
    Uy = -1.95*d(I,end).*n2(:,end);
    F(2*I(1:size(Y(I,:),1),1)-1,1) = Ux;
    F(2*I(1:size(Y(I,:),1),1),1)   = Uy;
    
    pDof = zeros(2*length(I),1);
    fDof = zeros(2*length(I),1);
    for kk = 1:length(I)
        pDof(2*kk-1) = 2*I(kk)-1; pDof(2*kk) = 2*I(kk);
        fDof(2*kk-1) = Ux(kk); fDof(2*kk) = Uy(kk);
    end
    Fem.PrescribedDisplacement = true;
end

if Fem.PrescribedDisplacement
    NLoad = size(Fem.Load,1);
    if ~isempty(Fem.Load)
        
        pDof = [];
        if numel(Fem.Load(1,2:end)) < 6
            for ii = 1:Fem.Dim
                pd = reshape(Fem.Dim*Fem.Load(:,1)+(ii-Fem.Dim),NLoad,1);
                F(pd) = Fem.Load(1:NLoad,1+ii);
                if norm(Fem.Load(1:NLoad,1+ii)) > 0
                    pDof = [pDof;pd];
                end
            end
        else
            R = reshape(Fem.Load(1,2:end),[Fem.Dim,Fem.Dim]);
            N0 = mean(Fem.Node(Fem.Load(:,1),:),1);
            V = (R*(Fem.Node(Fem.Load(:,1),:)-N0).').' - N0;
            for ii = 1:Fem.Dim
                pd = reshape(Fem.Dim*Fem.Load(:,1)+(ii-Fem.Dim),NLoad,1);
                F(pd) = V(:,ii);
                pDof = [pDof;pd];
            end
        end
    end
    
    EMod = Fem.Material.Emod();
   
    if Fem.Nonlinear
        KTtmp = full(Ktr);
        I = eye(Fem.Dim*Fem.NNode,Fem.Dim*Fem.NNode);
        KTtmp(pDof,:)    = 0*I(pDof,:);
        KTtmp(pDof,pDof) = -EMod*eye(length(pDof));
        Ktr = KTtmp;

        if isempty(Fem.Contact) || numel(Fem.Load(1,2:end)) >= 6
            alpha = Fem.TimeStep*Fem.OptFactor;
        else
            alpha = 1; 
        end
        
        if Fem.Iteration == 1
            F(pDof) = Fe(pDof)-EMod*F(pDof)*alpha;
        else 
            F(pDof) = Fe(pDof); 
        end
    else
        Ktmp = full(K); I = eye(Fem.Dim*Fem.NNode,Fem.Dim*Fem.NNode);
        Ktmp(pDof,:) = I(pDof,:);
        Ktmp(:,pDof) = I(:,pDof);
        F = -beta*K*F;
        if isempty(Fem.Contact), F(pDof) = beta*Fem.Load(1:NLoad,id+1);
        else, F(pDof) = beta*fDof;
        end
            
        K = Ktmp;
    end
end

if ~isempty(Fem.Contraction)
    Fb = sparse(Fem.i,Fem.e,Fem.fb0);
    F(:,1) = beta*Fb*Fem.Contraction(:);
end

if ~isempty(Fem.Pressure)
    IdMat = sign(Fem.Mesh.get('NodeToFace'));
    N = Fem.Normal;
    Elem = Fem.Element;
    PressureForce = zeros(length(Fem.Pressure),2);
    for ii = 1:length(Fem.Pressure)
        IdNode = Fem.Pressure(ii,1);
        Pressmag = max(Fem.Pressure(ii,2),Fem.Pressure(ii,3));
        el = find(IdMat(IdNode,:) > 1e-12);
        ntmp = [0,0];
        etmp = [];
        kk = 1;
        for jj = el
            w = Elem{jj}*0; w(Elem{jj} == IdNode) = 1;
            ntmp = w(:)'*N{jj};
            ee = Fem.Edge{jj};
            etmp(kk) = max(ee(Elem{jj} == IdNode));
            kk = kk+1;
        end
        
        ntmp = 0.5*(ntmp/length(el));
        PressureForce(ii,:) = ntmp*Pressmag;
    end
    
    NPLoad = size(Fem.Pressure,1);
    F(Fem.Dim*Fem.Pressure(1:NPLoad,1)-1,1) = beta*PressureForce(1:NPLoad,1)/NPLoad;
    F(Fem.Dim*Fem.Pressure(1:NPLoad,1),1)   = beta*PressureForce(1:NPLoad,2)/NPLoad;
end

if Fem.VolumetricPressure
    
    Pc = Fem.Mesh.get('Center');
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
    idnull = setdiff(1:Fem.NElem,id);
    W = ones(Fem.NElem,1);
    W(idnull) = 0;
    Ft = beta*sparse(Fem.i,1,W(Fem.e).*Fem.ft);
    
    if Fem.SolverStartMMA || Fem.Nonlinear
        z0 = Fem.Density;
        dz = Fem.ChangeMax;
        Fem.dFdE = zeros(2*Fem.NNode,Fem.NElem);
        bnd = boundary(Pc(id,1),Pc(id,2));
        idbound = id(bnd);
        for ii = idbound
            ze = z0;
            ze(ii) = ze(ii) + dz;
            Fem.Density = ze;
            id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
            idnull = setdiff(1:Fem.NElem,id);
            W = ones(Fem.NElem,1);
            W(idnull) = 0;
            fe = sparse(Fem.i,1,W(Fem.e).*Fem.ft);
            Fem.dFdE(:,ii) = ((fe - Ft)/(dz));
        end
        Fem.Density = z0;
    end
else
    Ft = sparse(Fem.NNode*Fem.Dim,1);
end

Fem.fInternal        = Fe; 
Fem.fExternal        = F - spMat*Fem.Utmp + Ft;
Fem.Stiffness        = K + spMat;
Fem.TangentStiffness = Ktr + spMat;
Fem.Residual         = Fem.fInternal - Fem.fExternal;
Fem.OutputVector     = L;

end

%---------------------------------------------- get free degrees-of-freedom
function FreeDofs = GetFreeDofs(Fem)
NSupp = size(Fem.Support,1);
if Fem.Dim == 2
FixedDofs = [Fem.Support(1:NSupp,2).*(2*Fem.Support(1:NSupp,1)-1);
             Fem.Support(1:NSupp,3).*(2*Fem.Support(1:NSupp,1))];
elseif Fem.Dim == 3
FixedDofs = [Fem.Support(1:NSupp,2).*(3*Fem.Support(1:NSupp,1)-2);
             Fem.Support(1:NSupp,3).*(3*Fem.Support(1:NSupp,1)-1);
             Fem.Support(1:NSupp,4).*(3*Fem.Support(1:NSupp,1))];
end
FixedDofs = unique(FixedDofs(FixedDofs>0));
AllDofs = 1:Fem.Dim*Fem.NNode;   
FreeDofs = setdiff(AllDofs,FixedDofs);
end

%--------------------------------------------------- local element matrices
function [Fe,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te] = Locals(Fem,eNode,eDof,dV,Rb)
% get order
nn = length(eNode);
mm = Fem.Dim;

% get gauss points and weights
W   = Fem.ShapeFnc{nn}.W;
Q   = Fem.ShapeFnc{nn}.Q;    
Fe  = zeros(mm*nn,1);
Te  = zeros(mm*nn,1);
Me  = zeros(mm*nn,mm*nn);
Ce  = zeros(mm*nn,mm*nn);
Ke  = zeros(mm*nn,mm*nn);
Kte = zeros(mm*nn,mm*nn);
SGP = zeros(length(W),6);
EGP = zeros(length(W),6);
Qe  = ones(nn,1);

Nshp = length(Fem.ShapeFnc{nn}.N(:,:,1));
NNe  = zeros(length(W),Nshp);

if mm == 2, Et = [dV;dV;0];
else, Et = [dV;dV;dV;0;0;0];
end

% quadrature loop
for q = 1:length(W)
    % extract shape-functions
    N     = Fem.ShapeFnc{nn}.N(:,:,q);
    dNdxi = Fem.ShapeFnc{nn}.dNdxi(:,:,q);
    J0    = Fem.Node0(eNode,:).'*dNdxi;
    dNdx  = dNdxi/J0;
    dJ    = det(J0);
       
    % get displacement field
    Delta = Fem.Utmp(eDof,:);
    
    % deformation gradient   
    F = DeformationGradient(Fem,Delta,dNdx);
    
    % polar decompostion
    [R,S,~] = PolarDecomposition(Fem,F);

    % increase robustness low density
    S = Rb*(S-eye(3)) + eye(3);
    
    % reconstruct deformation gradient
    F = R*S;
    
    % right cauchy-green strain
    C = F.'*F;
    
    if strcmp(Fem.Type,'PlaneStress') && Fem.Dim < 3
        C(3,3) = det(F)/(C(1,1)*C(2,2) - C(1,2)*C(2,1));
    end

    % get internal stress matrix
    [S0,D0] = Fem.Material.PiollaStress(C);
    
    % voigt-notation vectorize
    S = VoightNotation(S0);
    
    % reduced isotropic matrices
    [Se, De, Ge] = IsotropicReduction(Fem,D0,S);
    
    % nonlinear strain-displacement operator
    [Bnl,Bg,NN,tau] = NonlinearStrainOperator(Fem,N,dNdx,F);
    
    % internal force vector
    Fe = Fe + tau*W(q)*Bnl.'*Se*dJ;
    
    % lineararized stiffness matrix
    Ke = Ke + tau*W(q)*(Bnl.'*De*Bnl)*dJ;
    
    % tangent stiffness matrix
    Kte = Kte + tau*W(q)*(Bnl.'*De*Bnl + Bg.'*Ge*Bg)*dJ;
    
    % mass matrix
    Me = Me + tau*W(q)*Fem.Rho*(NN.')*NN*dJ;
    
    % dampings matrix
    Ce = Ce + tau*W(q)*Fem.Zeta*(NN.')*NN*dJ;
    
    % thermal expansion force
    Te = Te + tau*W(q)*Bnl.'*De*Et*dJ;
    
    % lagrangian strain
    Elagran = (1/2)*(C - eye(3));
    
    % true stress and lagrangian strain
    SGP(q,:) = VoightNotation((1/det(F))*F*S0*(F.'));
    EGP(q,:) = VoightNotation(Elagran);
    
    % construct shape functions
    NNe(((q-1)*Nshp + 1):(q*Nshp)) = N(:).';
end

SS = NNe.'*SGP;
EE = NNe.'*EGP;
[Svm, ~] = VonMises(SS(:,1),SS(:,2),SS(:,3),SS(:,4),SS(:,5),SS(:,6));
Svme = Svm(:); 

end

%----------------------------------------------- compute displacement field
function [ux,uy,un] = DisplacementField(Fem,U)
ux = zeros(Fem.NNode,1);
uy = zeros(Fem.NNode,1);
un = zeros(Fem.NNode,1);

for node = 1:Fem.NNode
    ux(node) = U(2*node - 1,1);
    uy(node) = U(2*node,1);
    un(node) = sqrt(ux(node)^2 + uy(node)^2);
end

end

%------------------------------------------------------ compute force field
function [fx,fy,fn] = InternalForceField(Fem,F)

fx = zeros(Fem.NNode,1);
fy = zeros(Fem.NNode,1);
fn = zeros(Fem.NNode,1);

for node = 1:Fem.NNode
    fx(node) = F(2*node - 1,1);
    fy(node) = F(2*node,1);
    fn(node) = sqrt(fx(node)^2 + fy(node)^2);
end

end

%----------------------------------------------------- deformation gradient
function F = DeformationGradient(Fem,U,dNdx)
nn = length(U)/Fem.Dim;
UU = zeros(nn,Fem.Dim);
UU(:,1) = U(1:Fem.Dim:Fem.Dim*nn);
UU(:,2) = U(2:Fem.Dim:Fem.Dim*nn);
if Fem.Dim == 2
F = (dNdx'*UU)';
F = [F(1,1)+1,F(1,2),0; F(2,1),F(2,2)+1,0;0,0,1];
else
UU(:,3) = U(3:Fem.Dim:Fem.Dim*nn);   
F = (dNdx'*UU)' + eye(3);
end
end

%---------------------------------------------------------- polar decompose
function [R,S,V] = PolarDecomposition(~,F)
[M, N] = size(F);
if M ~= N
    error('Matrix must be square.');
end
C = F'*F;
[Q0, lambdasquare] = eig(C);
lambda = sqrt(diag((lambdasquare))); 
Uinv = repmat(1./lambda',size(F,1),1).*Q0*Q0';
R = F*Uinv;
S = R'*F;
V = F*R';
end

%---------------------------------------------------------- strain operator
function [Bn,Bg,NN,tau] = NonlinearStrainOperator(Fem,N,dNdx,F)
nn = length(N);
mm = Fem.Dim;

NN = zeros(mm,mm*nn);
NN(1,1:mm:mm*nn) = N.';
NN(2,2:mm:mm*nn) = N.';
dNdxX = dNdx(:,1).';
dNdxY = dNdx(:,2).';

if Fem.Dim == 3
    NN(3,3:mm:mm*nn) = N.';
    dNdxZ = dNdx(:,3).';
end

Bn = zeros((mm-1)*3,mm*nn);
Bg = zeros((mm-1)*4+(mm-2),mm*nn);

if mm == 2 % 2-dimensional
    Bn(1,1:mm:mm*nn) = dNdxX*F(1,1);
    Bn(1,2:mm:mm*nn) = dNdxX*F(2,1);
    Bn(2,1:mm:mm*nn) = dNdxY*F(1,2);
    Bn(2,2:mm:mm*nn) = dNdxY*F(2,2);
    Bn(3,1:mm:mm*nn) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(3,2:mm:mm*nn) = dNdxX*F(2,2) + dNdxY*F(2,1);
    
    Bg(1,1:mm:mm*nn) = dNdxX;
    Bg(2,1:mm:mm*nn) = dNdxY;
    Bg(3,2:mm:mm*nn) = dNdxX;
    Bg(4,2:mm:mm*nn) = dNdxY;
    
else % 3-dimensional
    Bn(1,1:mm:mm*nn) = dNdxX*F(1,1);
    Bn(1,2:mm:mm*nn) = dNdxX*F(2,1);
    Bn(1,3:mm:mm*nn) = dNdxX*F(3,1);
    Bn(2,1:mm:mm*nn) = dNdxY*F(1,2);
    Bn(2,2:mm:mm*nn) = dNdxY*F(2,2);
    Bn(2,3:mm:mm*nn) = dNdxY*F(3,2);
    Bn(3,1:mm:mm*nn) = dNdxZ*F(1,3);
    Bn(3,2:mm:mm*nn) = dNdxZ*F(2,3);
    Bn(3,3:mm:mm*nn) = dNdxZ*F(3,3);
    Bn(4,1:mm:mm*nn) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(4,2:mm:mm*nn) = dNdxX*F(2,2) + dNdxY*F(2,1);
    Bn(4,3:mm:mm*nn) = dNdxX*F(3,2) + dNdxY*F(3,1);
    Bn(5,1:mm:mm*nn) = dNdxY*F(1,3) + dNdxZ*F(1,2);
    Bn(5,2:mm:mm*nn) = dNdxY*F(2,3) + dNdxZ*F(2,2);
    Bn(5,3:mm:mm*nn) = dNdxY*F(3,3) + dNdxZ*F(3,2);
    Bn(6,1:mm:mm*nn) = dNdxX*F(1,3) + dNdxZ*F(1,1);
    Bn(6,2:mm:mm*nn) = dNdxX*F(2,3) + dNdxZ*F(2,1);
    Bn(6,3:mm:mm*nn) = dNdxX*F(3,3) + dNdxZ*F(3,1);
    
    Bg(1,1:mm:mm*nn) = dNdxX;
    Bg(2,1:mm:mm*nn) = dNdxY;
    Bg(3,1:mm:mm*nn) = dNdxZ;
    Bg(4,2:mm:mm*nn) = dNdxX;
    Bg(5,2:mm:mm*nn) = dNdxY;    
    Bg(6,2:mm:mm*nn) = dNdxZ;   
    Bg(7,3:mm:mm*nn) = dNdxX;
    Bg(8,3:mm:mm*nn) = dNdxY;    
    Bg(9,3:mm:mm*nn) = dNdxZ;     
end

tau = 1;
end 

%------------------------------------------------------ isotropic reduction
function [S, D, G] = IsotropicReduction(Fem,D0,S0)
mm = Fem.Dim;
D = zeros((mm-1)*3,(mm-1)*3); 
S = zeros((mm-1)*3,1);

if mm == 2
    D(1:2,1:2) = D0(1:2,1:2); D(3,3) = D0(4,4);
    SIG = [S0(1),S0(4); S0(4),S0(2)];
    S(1) = S0(1); S(2) = S0(2);
    S(3) = S0(4);
else
    D = D0;
    S = S0;
    SIG = [S0(1),S0(4),S0(6); S0(4),S0(2),S0(5); S0(6), S0(5),S0(3)];
end

G = kron(eye(mm),SIG);
end

%------------------------------------------------------ isotropic reduction
function [flag,Fem] = CheckConvergence(Fem,SingularKt)

CriteriaResidual = (Fem.SolverResidual(end) > Fem.ResidualNorm);

DiffSVM = abs(Fem.SolverVonMises(end-1)/Fem.SolverVonMises(end) - 1);

CriteriaStress = (DiffSVM > Fem.StressNorm);

Criteria = CriteriaResidual && CriteriaStress;

if Fem.SolverResidual(end,1) > Fem.SolverResidual(end-1,1)
   Fem.Divergence = Fem.Divergence + 1;
end

if (Criteria && (Fem.Iteration <= Fem.MaxIteration) && ~SingularKt && ...
         Fem.Divergence < 5)
    flag = 0;
else
    if (Fem.Iteration > Fem.MaxIteration) || SingularKt || ...
            Fem.SolverResidual(end) > (25*Fem.Material.Emod)
        flag = 2; 
        Fem.Convergence = false;
        Fem.Utmp = Fem.U;
    elseif (Fem.Iteration == 1 && Fem.Nonlinear)
        flag = 0;    
    elseif (Fem.Divergence >= 5)
        flag = 2;        
        Fem.Utmp = Fem.U;
        Fem.Convergence = false;
    else
        flag = 1;
        Fem.Convergence = true;
    end
    if round((100*Fem.Iteration/Fem.MaxIteration)) > 75 && ~Fem.SolverStartMMA
        Fem.MaxIteration = round(Fem.MaxIteration)*1.25;
    end
end

ProcessMonitor(Fem);

end

%--------------------------------------------------- numerical solver timer
function [Fem, Terminate] = SolverTimer(Fem)
    
Terminate = false;
if ~Fem.Nonlinear, Fem.Time = 1-0.5*Fem.TimeStep;  end
    
if Fem.Convergence
    Fem.U = Fem.Utmp;
    if Fem.Time + Fem.TimeStep > Fem.TimeEnd  && ~Fem.EndIncrement ...
            && Fem.Time ~= Fem.TimeEnd
        Fem.TimeStep = Fem.TimeEnd-Fem.Time-1e-6;
        Fem.EndIncrement = true;
    end
    
    Fem.Time = Fem.Time + Fem.TimeStep;
    
    if Fem.Time > Fem.TimeEnd
        Terminate = true;
    end
    
elseif ~Fem.Convergence
    Fem.Utmp = Fem.U;
    Fem.TimeStep = clamp(Fem.TimeStep/2,Fem.TimeStepMin,1);
    Fem.Time = clamp(Fem.Time - Fem.TimeStep,0,1);
    Fem.BisectCounter = Fem.BisectCounter+1;
    fprintf('------------------------------------------------------------|\n');
    if ~Fem.SolverStartMMA
    Fem.MaxIteration = Fem.MaxIteration + 15;
    fprintf(' Bisection                                                  |\n');
    fprintf('------------------------------------------------------------|\n');
    end
    if Fem.BisectCounter > Fem.BisectLimit
        Terminate = true;
    end
end

if ~Fem.SolverStart, Fem.Increment = 1; Fem.SolverStart = true;
else, Fem.Increment = Fem.Increment + 1;
end

if (Fem.SolverStartMMA && Terminate)||(Fem.SolverStartMMA && ~Fem.Nonlinear)
    %ProcessMonitor(Fem); 
end

end

%----------------------------------------------------------- mesh smoothing
function f = Smoothing(Fem,f,naver)
face = Fem.Element; 
W = NodeAdjecency(face);
n = length(W); 
D = spdiags(full(sum(W,2).^(-1)),0,n,n);
W = D*W;

for ii=1:naver
    f = W*f(:);
end

end

%%/////////////////////////////////////////////////// TOPOLOGY OPTIMIZATION
%------------------------------------------------------ objective function
function [f,dfdE,dfdV] = ObjectiveFunction(Fem)
    
u(:,1) = Fem.Utmp;
fDof = GetFreeDofs(Fem);
E = MaterialField(Fem);

if strcmp(Fem.OptimizationProblem,'Compliance') &&  ~Fem.Nonlinear
f = Fem.fExternal.'*u;
temp = cumsum(-u(Fem.i).*Fem.k.*u(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliance') &&  Fem.Nonlinear
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
f = Fem.fExternal.'*u;
lam = 0*u(:,1);
lam(fDof) = K(fDof,fDof)\Fem.fExternal(fDof);
temp = cumsum(-u(Fem.i).*Fem.k.*lam(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliant') && ~Fem.Nonlinear
E = MaterialField(Fem);
f = -Fem.OutputVector.'*u(:,1); 
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
lam = 0*u(:,1);
lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
temp = cumsum(u(Fem.i,1).*Fem.k.*lam(Fem.j,1));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliant') && Fem.Nonlinear
E = MaterialField(Fem);
f = -Fem.OutputVector.'*u(:,1); 
%K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
lam = 0*u(:,1);
lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
temp = cumsum(u(Fem.i).*Fem.k.*lam(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
end

if Fem.VolumetricPressure 
temp = cumsum(Fem.dFdE(Fem.e).*lam(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = dfdE + [temp(1);temp(2:end)-temp(1:end-1)];
end

dfdV = zeros(size(Fem.NElem));

end

%------------------------------------------------------ constraint function
function [g,dgdE,dgdV] = ConstraintFunction(Fem)
V = Fem.SpatialFilter*Fem.Density;
A = Fem.Mesh.get('Area');
g = sum(A.*V)/(sum(A)*Fem.VolumeInfill)-1;
dgdE = zeros(size(A));
dgdV = A/(sum(A)*Fem.VolumeInfill);
end

%---------------------------------------------- method of moving asymptotes
function [Fem,zNew] = UpdateSchemeMMA(Fem,f,dfdz,g,dgdz)  
iter = Fem.IterationMMA; 
N = length(Fem.Density);
M = 1;

if isempty(Fem.fnorm), Fem.fnorm = abs(norm(dfdz)); end

xval = Fem.Density;
xmin = zeros(N,1);
xmax = ones(N,1);
f0val = (f/norm(f));
df0dx = (dfdz/norm(dfdz));
df0dx2 = 0*df0dx;
fval = g;
dfdx = dgdz;
dfdx2 = dgdz*0;

A0 = 1; A = 0; C = 10000*ones(M,1); D = 0;

if iter == 1, Fem.low = xmin; Fem.upp = xmax; end
if iter > 1, Fem.xold1 = Fem.xold1; else, Fem.xold1 = 0; end
if iter > 2, Fem.xold2 = Fem.xold2; else, Fem.xold2 = 0; end

[xmma,~,~,~,~,~,~,~,~,Fem.low,Fem.upp] = mmasub(M,N,iter,xval,xmin,xmax,...
    Fem.xold1,Fem.xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,Fem.low,...
    Fem.upp,A0,A,C,D);

if strcmp(Fem.OptimizationProblem,'Compliance')
dx = clamp(xmma - xval,-Fem.ChangeMax*15,Fem.ChangeMax*15);
else
dx = clamp(xmma - xval,-Fem.ChangeMax*15,Fem.ChangeMax*15);    
end
zNew = xval + dx;

if iter >= 1, Fem.xold1 = xval; end
if iter >= 2, Fem.xold2 = Fem.xold1; end

end

%---------------------------------------------- method of moving asymptotes
function [Fem,zNew] = UpdateSchemeOC(Fem,dfdz,g,dgdz)  
move=Fem.ChangeMax;
eta=0.5;
l1=0; l2=1e6;  
z0 = Fem.Density;

if strcmp(Fem.OptimizationProblem,'Compliance')
    while l2-l1 > 1e-4
        lmid = 0.5*(l1+l2);
        B = -(dfdz./dgdz)/lmid;
        zCnd = Fem.zMin+(z0-Fem.zMin).*B.^eta;
        zNew = max(max(min(min(zCnd,z0+move),Fem.zMax),z0-move),Fem.zMin);
        if (g+dgdz'*(zNew-z0)>0),l1=lmid;
        else, l2=lmid; end
    end
    
elseif strcmp(Fem.OptimizationProblem,'Compliant')
    while (l2-l1)/(l2+l1) > 1e-4 && l2>1e-40
        lmid = 0.5*(l1+l2);
        B = max(1e-10,-(dfdz./dgdz)/lmid);
        zCnd = Fem.zMin+(z0-Fem.zMin).*B.^eta;
        zNew = max(max(min(min(zCnd,z0+move),Fem.zMax),z0-move),Fem.zMin);
        if (g+dgdz'*(zNew-z0)>0),l1=lmid;
        else, l2=lmid;   end
    end
end
end

%------------------------------------------------------ isotropic reduction
function [bool,Fem] = CheckConvergenceOpt(Fem)

if ((norm(1) > 1e-3) && (Fem.IterationMMA <= Fem.MaxIterationMMA))
    bool = true;
    if (mod(Fem.IterationMMA,Fem.PenalStep) == 0 && ...
            Fem.Penal < Fem.PenalMax)
        Fem.Penal = Fem.Penal + 1;
    end
else
    bool = false;
end

end

%--------------------------------------------------- material interpolation 
function [E,dEdy,V,dVdy] = MaterialField(Fem,ForceFilterOff)
if nargin < 2, ForceFilterOff = false; end

if ~ForceFilterOff, y = Fem.SpatialFilter*Fem.Density;
else, y = Fem.Density; end

eps = Fem.Ersatz;
switch(Fem.MaterialInterpolation)
  case('SIMP')
    E = eps+(1-eps)*y.^Fem.Penal;
    V = y;
    dEdy = (1-eps)*Fem.Penal*y.^(Fem.Penal-1);
    dVdy = ones(size(y,1),1);
  case('SIMP-H')
    penal = Fem.Penal;
    beta = Fem.Beta; 
    h = 1-exp(-beta*y)+y*exp(-beta);
    E = eps+(1-eps)*h.^penal;
    V = h;
    dhdy = beta*exp(-beta*y)+exp(-beta);
    dEdy = (1-eps)*penal*h.^(penal-1).*dhdy;
    dVdy = dhdy;
end
end

%---------------------------------------------------------- material filter 
function P = GenerateRadialFilter(Fem,R)
if R<0, P = speye(Fem.NElem); return; end

PS1 = Fem.Mesh.get('Center');
PS2 = Fem.Mesh.get('Center');

d = DistancePointSet(Fem,PS1,PS2,R);

P = sparse(d(:,1),d(:,2),1-d(:,3)/R);
P = spdiags(1./sum(P,2),0,Fem.NElem,Fem.NElem)*P;
end

%---------------------------------------------------------- material filter 
function Z = InitialDesign(Fem,Arg)

N = Arg{1}; Nx = N(1); Ny = N(2);
R = Arg{2};

Pc = Fem.Mesh.get('Center');

dX = (Fem.BdBox(2)-Fem.BdBox(1)-2*R*Nx)/(1+Nx);
dY = (Fem.BdBox(4)-Fem.BdBox(3)-2*R*Ny)/(1+Ny);

dXX = 0:(Nx-1); dXX = (dX + R)+dXX*(2*R + dX);
dYY = 0:(Ny-1); dYY = (dY + R)+dYY*(2*R + dY);

[Xc,Yc] = meshgrid(dXX,dYY);

CenCircle = [Xc(:),Yc(:)];

d = DistancePointSet(Fem,CenCircle,Pc,R);

Z = ones(Fem.NElem,1); 
Z(d(:,1)) = 0;

end

%---------------------------------------------------------- material filter 
function d = DistancePointSet(Fem,PS1,PS2,R)
B = Fem.BdBox;
d = cell(size(PS1,1),1);
B1 = B(1); B2 = B(2); B12 = lerp(B1,B2,0.5);
B3 = B(3); B4 = B(4); B34 = lerp(B3,B4,0.5);
if Fem.Dim == 3
B5 = B(5); B6 = B(6); B56 = lerp(B5,B6,0.5);
end
    
for el = 1:size(PS1,1)   
    if Fem.Dim == 3
    dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2 + ...
       (PS1(el,3)-PS2(:,3)).^2 );    
    else    
    dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    end
    
    if ~isempty(Fem.Periodic)
        Rp = Fem.Periodic;
        if Rp(1) == 1
             dist2 = sqrt((PS1(el,1)-(PS2(:,1) + B2-B1)).^2 +...
             (PS1(el,2)-PS2(:,2)).^2);
             dist = (min([dist,dist2].'))';
        end
        if Rp(1) == 1/2
            PS = -(PS2(:,1)-B12) + B12;
            dist2 = sqrt((PS1(el,1)-PS).^2 +...
             (PS1(el,2)-PS2(:,2)).^2);
             dist = (min([dist,dist2].'))';
        end
        if Rp(2) == 1/2
            PS = -(PS2(:,2)-B34) + B34;
            dist2 = sqrt((PS1(el,1)-PS2(:,1)).^2 +...
             (PS1(el,2)-PS).^2);
             dist = (min([dist,dist2].'))';
        end
    end

    [I,J] = find(dist<=R);   
    d{el} = [I,J+(el-1),dist(I)];
end

% matrix of indices and distance value
d = cell2mat(d); 
end

%--------------------------------------------------------------------------
function Rp = Reflection(Fem,P)
Domain = @(x) dRectangle(x,Fem.BdBox(1),Fem.BdBox(2),...
    Fem.BdBox(3),Fem.BdBox(4));
Area = (Fem.BdBox(2)-Fem.BdBox(1))*(Fem.BdBox(4)-Fem.BdBox(3));
d = Domain(P);

eps = 1e-8; 
eta = 0.95;

Alpha = 4*sqrt(Area/length(d));

NBdrySegs = size(d,2)-1;       
n1 = (Domain(P+repmat([eps,0],length(d),1))-d)/eps;
n2 = (Domain(P+repmat([0,eps],length(d),1))-d)/eps;

I  = abs(d(:,1:NBdrySegs))<Alpha; 

P1 = repmat(P(:,1),1,NBdrySegs);  
P2 = repmat(P(:,2),1,NBdrySegs); 
Rp(:,1) = P1(I)-2*n1(I).*d(I);  
Rp(:,2) = P2(I)-2*n2(I).*d(I);

d_R_P = Domain(Rp);
J    = abs(d_R_P(:,end))>=eta*abs(d(I)) & d_R_P(:,end)>0;
Rp   = Rp(J,:); 
Rp   = unique(Rp,'rows');
end

%-------------------------------------------------------------- GEN MESH
function SDF = GenerateCell(Fem,SDF0)
vq = SDF0;
    
if ~isempty(Fem.ReflectionPlane)
    RP = Fem.ReflectionPlane;
    if RP(1) == 1 && RP(2) == 1
        V = flip(vq,2); V = horzcat(V,vq);
        VV = flip(V,1); SDF = vertcat(VV,V);
    elseif RP(1) == 1 && RP(2) ~= 1
        V = flip(vq,2); V = horzcat(vq,V);
        SDF = V;
    elseif RP(1) == -1 && RP(2) ~= 1
        V = flip(vq,2); V = horzcat(V,vq);
        SDF = flip(V);
    elseif RP(1) ~= 1 && RP(2) == 1
        V = flip(vq,1); SDF = vertcat(V,vq);
    elseif RP(1) ~= 1 && RP(2) == -1
        V = flip(vq,1); SDF = vertcat(vq,V);
    elseif RP(1) == 0 && RP(2) == 0
        SDF = flip(vq);
    end
else
    SDF = flip(vq);
end

end

end
end

%---------------------------------------------------------- material filter 
function ProcessMonitor(Fem)
FreeDofs = GetFreeDofs(Fem);

if Fem.SolverStartMMA
if (Fem.IterationMMA == 1 && Fem.Iteration == 1 && Fem.Increment == 1)
fprintf(' Itr  | Inc  | k     | Obj. fun. | g(x)  | Residual  | Delta |\n');
fprintf('--------------------------------------------------------------\n');
end

fprintf(' %1.0f\t | %1.0f\t | %1.0f\t | %1.3e | %0.3f | %1.3e | %0.3f |\n',...
  Fem.IterationMMA,Fem.Increment,Fem.Iteration,abs(Fem.Objective(end)),...
  Fem.Constraint(end)+1,norm(Fem.Residual(FreeDofs)),norm(Fem.Change));

else
if Fem.Iteration == 1 && Fem.Increment == 1
fprintf(' Inc | Iter  | Residual  | Max. Svm  | Time | Beta  | dt     |\n');
fprintf('--------------------------------------------------------------\n');
end

fprintf(' %1.0f\t | %1.0f\t | %1.3e | %1.3e | %1.2f | %1.3f | %1.3f  |\n',...
    Fem.Increment,Fem.Iteration,norm(Fem.Residual(FreeDofs)),...
    max(Fem.s(:,1)),Fem.Time,Fem.LoadingFactor,Fem.TimeStep);   

end 

end

%------------------------------------------------------ optimization solver
function showInformation(Fem,Request)
if nargin < 2, Request = ''; end
    
if ~Fem.InformationBoolean
NodeNum = Fem.NElem;
MaxNVer = max(cellfun(@numel,Fem.Element));      
MinNVer = min(cellfun(@numel,Fem.Element));      

switch(Request)
case('NonlinearFem')
fprintf('--------------------------------------------------------------\n');  
if Fem.NElem>1e3, fprintf('* Elements: %1.2fk \n',NodeNum/1e3);
else, fprintf('* Element = %i \n',NodeNum); end
fprintf('* Element degree = P%1.0f-P%1.0f \n',MinNVer,MaxNVer);
fprintf('* Max iteration = %i \n', Fem.MaxIteration);
if Fem.Nonlinear, fprintf('* Nonlinear geometric = true\n');
else, fprintf('* Nonlinear geometric = false \n'); end
showMaterialInfo(Fem)
fprintf('--------------------------------------------------------------\n');
    
case('TopologyOptimization')
fprintf('--------------------------------------------------------------\n');  
if Fem.NElem>1e3, fprintf('* Elements: %1.2fk \n',NodeNum/1e3);
else, fprintf('* Element = %i \n',NodeNum); end
fprintf('* Element degree = P%1.0f-P%1.0f \n',MinNVer,MaxNVer);
fprintf('* Filter radius = %1.2f \n', Fem.FilterRadius);
fprintf('* Max iteration = %i \n', Fem.MaxIterationMMA);
fprintf('* Optimization = %s \n', Fem.OptimizationProblem);
fprintf('* Material scheme = %s \n', Fem.MaterialInterpolation);
if Fem.Nonlinear, fprintf('* Nonlinear geometric = true\n');
else, fprintf('* Nonlinear geometric = false\n'); end
showMaterialInfo(Fem)
fprintf('--------------------------------------------------------------\n');
end

Fem.InformationBoolean = true;
end

end

%------------------------------------------------- show material properties
function showMaterialInfo(Fem)
fprintf('--------------------------------------------------------------\n');
if strcmp(Fem.Material.Type,'Mooney')
fprintf('* Material model = Mooney-Rivlin (2-nd order) \n');
c1 = Fem.Material.C10; c2 = Fem.Material.C01; d = Fem.Material.K;
fprintf('\tC10 = %1.e, C01 = %1.3e, K = %1.3e \n', c1,c2,d);
elseif strcmp(Fem.Material.Type,'Yeoh')
fprintf('* Material model = Yeoh (3-nd order) \n');
c1 = Fem.Material.C1*1e3; c2 = Fem.Material.C2*1e3; c3 = Fem.Material.C3*1e3;
fprintf('\tC1 = %1.1e kPa, C2 = %1.1e kPa, C3 = %1.1e kPa\n', c1,c2,c3);
elseif strcmp(Fem.Material.Type,'NeoHookean')
fprintf('* Material model = Neo-Hookean \n');
E = Fem.Material.E*1e3; Nu = Fem.Material.Nu; 
fprintf('\tE = %1.1e kPa, Nu = %1.2f [-] \n', E,Nu);
end
end

%-------------------------------------------------------------- movie maker
function MovieMaker(Fem,Name,Request)
if nargin < 2, Request = ''; end

if Fem.Movie
    switch(Request)
        case('Start')
            filename = string([Name,'_', char(datetime(now,...
                'ConvertFrom','datenum')),'.gif']);
            
            filename = erase(filename,[":"," "]);
            background(gitpage);
            if ~isempty(Fem.MovieAxis), axis(Fem.MovieAxis); end
            if ~isempty(Fem.MovieCAxis), caxis(Fem.MovieCAxis); end
            drawnow;
            gif(char(filename),'frame',gcf,'nodither');
        otherwise
            background(gitpage);
            if ~isempty(Fem.MovieAxis), axis(Fem.MovieAxis); end
            if ~isempty(Fem.MovieCAxis), caxis(Fem.MovieCAxis); end
            drawnow;
            gif;
    end
end

end

%----------------------------------------------- triangulate polygonal mesh
function [v, f] = QuickTriangulation(Fem,Center,v0,f0)
f = [];

for ii = 1:Fem.NElem
    el = f0{ii};
    n = numel(el);
    elem = [el(1:n)', [el(2:n)'; el(1)], ...
            repmat(Fem.NNode+ii,n,1)];
    f = [f; elem];
end

f = num2cell(f,2);
v = [v0;Center];
end

%------------------------------------------------ generate elemental matrix
function ElemMat = GenerateElementalMatrix(Element,NElem)
El = Element(1:NElem)';                 
MaxNVer = max(cellfun(@numel,El));      
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))]; 
ElemMat = cellfun(PadWNaN,El,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});       
end
