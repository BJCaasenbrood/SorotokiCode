% class Fem(msh,varargin)
%--------------------------------------------------------------------------
% FEM is a class used for (nonlinear) finite element methods and topology
% optimization.  
% -------------------------------------------------------------------------  
% main usage:
%   fem = Fem(msh);                      % converts MESH class to FEM class
% ------------------------------------------------------------------------- 
% options:
%   fem = Fem(msh,'Timestep',1e-3);                   % set time increments
%
% -------------------------------------------------------------------------  
% add boundary constraints:
%   fem = fem.AddConstraint('Support',id,[1,0]);     % fixes x-displacement
%   fem = fem.AddConstraint('Support',id,[0,1]);     % fixes y-displacement 
%   fem = fem.AddConstraint('Support',id,[1,1]);     % no displacement x-y
%
%   fem = fem.AddConstraint('Load',id,[A,0]);      % apply load A in x-axis
%   fem = fem.AddConstraint('Load',id,[0,A]);      % apply load A in y-axis
%   (Note: id are Node indices)
%
%   fem = fem.AddConstraint('Displace',id,[A,0]);    % displace A in x-axis
%   fem = fem.AddConstraint('Displace',id,[0,A]);    % displace A in y-axis
%   (Note: id are Node indices)
%
%   fem = fem.AddConstraint('Pressure',id,A);            % apply pressure A 
%   (Note: id are Edge/Surface indices)
% -------------------------------------------------------------------------  
% assign material:
%   fem.Material = YeohMaterial('C1',1,'C2',0.1);
%   fem.Material = Ecoflex0030;
%   fem.Material = Dragonskin10;
% 
%  Also see: SOROTOKI, MESH, SDF
% -------------------------------------------------------------------------    
classdef Fem < handle
    
    properties (Access = public)
        BdBox;        % Bounding box of FEM
        Mesh;         % Mesh class
        Material;     % Material class
        Density;      % Density field (Topology)
        Node;         % Vertices/Nodes
        Element;      % Elements
        NNode;        % Number of nodes
        NElem;        % Number of elements
        Dim;          % Dimension (2/3)
        Center;       % Element centers
        Topology;     % current save of topology
        Log;          % Log structure
        Time;         % Time (solver)
        Flow;         % auxilary flow function
        Contraction;  % list of elemental contractions
        
        Utmp = 0;
    end
    
    properties (Access = private)
        Node0;         % list of nodes in reference configuration
        Center0;       % list of element centers in reference configuration
        Support;       % list of boundary conditions (i.e., supports)
        Load;          % list of applied forces
        Tendon;        % list of applied follower forces (i.e., tendons)
        Spring;        % list of sping attachments to nodes
        Output;        % list of nodes that are logged
        Pressure;      % list of edges subjected to pressures
        PressureCell;  % list of elements undergoing volumetric expansion
        Contact;       % SDF contact function
        
        ElemNDof;      % Degrees-of-freedom of elements
        ElemSort;      % Elements in sorted list
        ElemRot;       % Rotation mat. each element
        ElemStr;       % Stretch vector each element
        Normal;        % Normals of each elements
        Edge;          % Cell of edges/boundaries
        Rotation;      % F = RQ: rotational part
        Stretch;       % F = RQ: stretch part
        
        PrescribedDisplacement = false;
        VolumetricPressure     = false;
        TendonForces           = false;
        GravityRamp            = true;
        
        VonMisesNodal; 
        sxxNodal; syyNodal; sxyNodal;
        exxNodal; eyyNodal; exyNodal;
        fxNodal;  fyNodal;
        rotNodal;
        
        SPxyNodal;
        Potential;
        PotentialG;
        PotentialF;
        Kinetic;
        ElemMat;
        Gravity;
        
        fInternal = rand(1)*1e-8;
        fExternal = rand(1)*1e-8;
        fInput    = rand(1)*1e-8;
        
        Stiffness; 
        TangentStiffness;
        MassMatrix;
        DampingMatrix;
        
        Residual = Inf;
        
        TimeEnd      = 1.0;
        TimeStep     = 0.1;
        TimeStep0    = 0.1;
        TimeDelta    = 0;
        TimeStepMin  = 1e-3;
        EndIncrement;
        LoadingFactor;
        SigmoidFactor;
        ShowProcess;
        
        ResidualNorm = 1e-5;
        StressNorm   = 1e-9;
        DisplaceNorm = 1e-9;
        Objective    = 1e-4;
        Constraint   = 1e-4;
        
        Type     = 'PlaneStrain'
        Solver   = 'mr';
        SolverId;
        
        IterationMMA = 1;
        Iteration    = 1;
        Increment    = 1;
        Divergence   = 0;
        Filter       = (1/5); 
        ShapeFnc;
        
        U        = 0;
        
        dUtmp    = 0;
        ddUtmp   = 0;
        Ptmp     = 0;
        dPtmp    = 0;
        ddPtmp   = 0;
        Ftmp     = 0;
        Utmp_    = 0;
        Ptmp_    = 0;
        dUtmp_   = 0;
        ddUtmp_  = 0;
        Z0       = 0;
        
        PlTmp = 0;
        
        MaxIteration     = 500;
        MaxIterationMMA  = 50;
        SolverStart      = false;
        SolverStartMMA   = false;
        SolverPlot       = true;
        SolverPlotType   = 'Svm';
        Assemble         = false;
        Convergence      = false; 
        Nonlinear        = true;
        BisectLimit      = 5;
        BisectCounter    = 1;
        AssembledSystem  = false;
        PressureLoad     = 0;
        
        Linestyle  = '-';
        Linestyle0 = '-';
        Colormap    = turbo;
        ColormapOpt = barney(-1);
        ColorAxis   = [];
        I3 = eye(3); 
        O3 = zeros(3);
        i; j; m; fi; t; e; c; s; p; v; l; k; fb; ed; fb0; ft;
        
        SolverResidual = 1e7;
        SolverVonMises = 1e7;
        SolverDisplace = 1e7;
        
        InformationBoolean = false;
        
        Movie = false;
        MovieStart = false;
        MovieAxis; MovieCAxis;
        
        VolumeInfill = 0.3;
        Ersatz       = 1e-3; 
        Penal        = 1;
        PenalMax     = 4;
        PenalStep    = 10;
        ChangeMax    = 0.03;
        Beta         = 1.5;
        FilterRadius = 0.5;
        SpatialFilter;
        
        OptFactor             = 1;
        OptimizationSolver    = 'MMA';
        MaterialInterpolation = 'SIMP';
        OptimizationProblem   = 'Compliance';
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
    obj.Center   = Mesh.get('Center');
    obj.Center0  = Mesh.get('Center');
    obj.NNode    = Mesh.get('NNode');
    obj.Dim      = Mesh.get('Dim');
    obj.Element  = Mesh.get('Element');
    obj.NElem    = Mesh.get('NElem');
    obj.Center   = Mesh.get('Center');
    obj.BdBox    = Mesh.get('BdBox');
    obj.Density  = ones(obj.NElem,1);
    obj.Residual = zeros(obj.Dim*obj.NNode,1);
    obj.Utmp     = zeros(obj.Dim*obj.NNode,1);
    obj.dUtmp    = zeros(obj.Dim*obj.NNode,1);
    obj.Gravity  = zeros(obj.Dim,1);
    obj.Flow     = [];
    
    obj.SigmoidFactor  = 0;
    obj.ShowProcess    = true;
    
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
        if strcmp(varargin{ii},'TimeStep') % overwrite initial t-step
            Fem.TimeStep0 = varargin{ii+1}; 
        end
        if strcmp(varargin{ii},'FilterRadius') % regenerate filter
            Fem.SpatialFilter = GenerateRadialFilter(Fem,varargin{ii+1});
        end
    end
    
end
%--------------------------------------------------------------------- show
function h = show(Fem,varargin)
% ensure a show selection    
if nargin<2
    Request = '0'; 
else
    Request = varargin{1}; 
end

flag = 0;
S    = 'interp'; 
V    = Fem.Node;
colormap(Fem.Colormap);

switch(Request)
    case('0'),     Z = zeros(Fem.NNode,1);
    case('SDF'),   Z = Fem.Mesh.SDF(Fem.Mesh.Node); Z = Z(:,end);
    case('Svm'),   Z = Fem.VonMisesNodal;
    case('Sxx'),   Z = Fem.sxxNodal;
    case('Syy'),   Z = Fem.syyNodal;
    case('Sxy'),   Z = Fem.sxyNodal;
    case('Exx'),   Z = Fem.exxNodal;
    case('Eyy'),   Z = Fem.eyyNodal;
    case('Exy'),   Z = Fem.exyNodal;
    case('Fx'),    Z = Fem.fxNodal;
    case('Fy'),    Z = Fem.fyNodal;
    case('Rot'),   Z = Fem.rotNodal;
    case('Field'), Z = varargin{2}; varargin = varargin{1};
    case('Fin'),   [~,~,Z] = DisplacementField(Fem,Fem.fInternal);
    case('Fex'),   [~,~,Z] = DisplacementField(Fem,Fem.fExternal);
    case('Un'),    [~,~,Z] = DisplacementField(Fem,Fem.Utmp);
    case('Ux'),    [Z,~,~] = DisplacementField(Fem,Fem.Utmp);
    case('Uy'),    [~,Z,~] = DisplacementField(Fem,Fem.Utmp);
    case('E'),     [~,~,Z] = MaterialField(Fem); S = 'flat'; 
                   V = Fem.Node0; colormap(barney(-1)); background('w');
    case('E+'),    [~,~,Z] = MaterialField(Fem); S = 'flat'; 
                   colormap(noir(-1)); background('w');
    otherwise;   flag = 1; Z = 0;
end

% ensures smooth plot
if length(Z) == Fem.NNode
    Z = Smoothing(Fem,Z,1); 
end

if length(varargin) == 2 && flag == 0
    handle = varargin{2}; 
    set(handle{2},'FaceVertexCData',Z); 
    set(handle{2},'Vertices',V); 
    return; 
end

if flag == 0
cla; axis equal;     
axis off; hold on; h{3} = [];

if ~isempty(Fem.ColorAxis)
   caxis(Fem.ColorAxis); 
   colorbar;
end

FaceMatrix  = Fem.Mesh.get('ElemMat');
BoundMatrix = Fem.Mesh.get('Boundary');

if Fem.Dim > 2 && (strcmp(Request,'E') || strcmp(Request,'E+'))
    Z(Z<Fem.VoidTolerance,:) = nan;
    S = 'flat';
    Z = Fem.Mesh.get('ElementToFace')*Z;
else
    Alp = ones(Fem.NElem,1);
end

h{1} = patch('Faces',FaceMatrix,'Vertices',Fem.Node0,...
    'LineStyle','none','Linewidth',1.5,'FaceColor',[0.9,0.9,0.9]);

h{2} = patch('Faces',FaceMatrix,'Vertices',V,...
    'FaceVertexCData',Z,'Facecolor',S,'LineStyle',Fem.Linestyle,...
    'Linewidth',1.5,'FaceAlpha',1,'EdgeColor','k');

h{3} = patch('Faces',BoundMatrix,'Vertices',V,...
    'LineStyle','-','Linewidth',2,'EdgeColor','k');

if Fem.Dim == 3, view(30,10); end

if ~isempty(Fem.Contact) && Fem.Dim == 2
    sdfFNC = Fem.Contact{1};
    Move = Fem.Contact{2};
    BD = Fem.BdBox;
    %BD = [-BD,BD,-BD,BD];
    [px,py] = meshgrid(linspace(BD(1),BD(2),50),linspace(BD(3),BD(4),50));
    Y = [px(:),py(:)];
    beta = 0.975*Fem.Time;
    Y(:,1) = Y(:,1) - beta*Move(1);
    Y(:,2) = Y(:,2) - beta*Move(2);
    d = sdfFNC(Y);
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
    case('ISO')
        clf; former(Fem);
        if nargin < 3, ISOVALUE = 0.2; else, ISOVALUE = varargin{2};
        end
        showISO(Fem,ISOVALUE,0.5);
    end
end

if Fem.Movie 
    background(metropolis);
    if Fem.MovieStart == false
       Fem.MovieStart = true;
       
       if ~Fem.SolverStartMMA
            Name = 'fem'; 
       else 
            Name = 'topo'; 
       end
       
       MovieMaker(Fem,Name,'Start');
    else
       MovieMaker(Fem);
    end
end

end
%-------------------------------------------------------------------- reset
function Fem = reset(Fem,varargin)
% reset simulation data,log    
if isempty(varargin), Request = 'fem';
else, Request = varargin{1};
end

switch(Request)
    case('fem')
        Fem.Iteration = 1;
        Fem.Increment = 1;
        %Fem.Node      = Fem.Node0;
        Fem.Support   = [];
        Fem.Load      = [];
        Fem.Spring    = [];
        Fem.Log       = [];
        Fem.Residual  = zeros(Fem.Dim*Fem.NNode,1);
        Fem.Utmp      = zeros(Fem.Dim*Fem.NNode,1);
        Fem.dUtmp     = zeros(Fem.Dim*Fem.NNode,1);
    otherwise
        Fem.Iteration    = 1;
        Fem.Increment    = 1;
        Fem.IterationMMA = 1;
        Fem.Node         = Fem.Node0;
        Fem.Utmp         = [];
        Fem.Change       = [];
        Fem.Objective    = [];
        Fem.Constraint   = [];
        Fem.Log          = [];
end

Fem = SetupFiniteElement(Fem);

end
%-------------------------------------------------- NL newton-raphson solve
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


TempNode          = cell(1);
TempNode{1}       = Fem.Node;

if Fem.ShowProcess
    showInformation(Fem,'NonlinearFem');
end

if Fem.SolverPlot || (Fem.SolverStartMMA && Fem.SolverPlot)
    figure(101); 
    Fem.show('0');
    clf;
    background(metropolis);
end

if norm(Fem.Utmp) == 0
    Fem.Utmp  = zeros(Fem.Dim*Fem.NNode,1);
    Fem.Node  = Fem.Node0;
else
    % update mesh class
    Fem.GravityRamp = false;
    Fem.Node    = UpdateNode(Fem,Fem.Utmp);
    Fem.Center  = UpdateCenter(Fem);
end

while true 
    
    [Fem, Terminate] = SolverTimer(Fem);
    if Terminate, break; end
        
    flag           = 0;
    Singular       = false;
    Fem.Divergence = 0;
        
    % load increment
    Fem.LoadingFactor = sigmoid(Fem.Time,Fem.SigmoidFactor);
    
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
            b = Fem.Residual(FreeDofs);
        elseif ~Fem.Nonlinear
            A = Fem.Stiffness(FreeDofs,FreeDofs);
            b = sparse(Fem.fExternal(FreeDofs) + Fem.fInput(FreeDofs));
        end
        
        Delta = Fem.Utmp;

        %if rcond(full(A)) >= 1e-10
            if Fem.SolverStartMMA
                %[L,D,P] = ldl(A,'vector');
                DeltaU = A\b; % topology optimization needs exact mid-solution
                %DeltaU = sparse(P,1,(L'\(D\(L\b(P))))); % thanks Ondrej ;)
            else
                
              % [L,D,P] = ldl(A,'vector');
               %minDiag = full(min(diag(D)));
               
%                if Fem.SolverId == 1
%                     DeltaU = A\b;
%                     %DeltaU = sparse(P,1,(L'\(D\(L\b(P))))); % thanks Ondrej ;)
%                elseif Fem.SolverId == 2
                    %DeltaU = sparse(P,1,(L'\(D\(L\b(P))))); % thanks Ondrej ;)
%                elseif Fem.SolverId == 3
%                    [DeltaU,~] = gmres(A,b,[],Fem.ResidualNorm,100);
%                else
                    DeltaU = A\b;
%                end

            end     
            
        if Fem.Nonlinear
            Delta(FreeDofs,1) = Delta(FreeDofs,1)-DeltaU(:,1);
        else
            Delta(FreeDofs,1) = DeltaU(:,1); 
            b = Fem.ResidualNorm; 
        end
        
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
            Fem.Utmp   = Delta;
            Fem.Node   = UpdateNode(Fem,Delta);
            Fem.Center = UpdateCenter(Fem);
        elseif flag == 1 && Fem.Nonlinear
            Fem.Utmp   = Delta;
            Fem.Node   = UpdateNode(Fem,Delta);
            Fem.Center = UpdateCenter(Fem);
        end
        
        % update mesh class
        Fem.Mesh = Fem.Mesh.set('Node',Fem.Node);
        Fem.Mesh = Fem.Mesh.set('Center',Fem.Center);

    end 
    
    if ~Fem.SolverStartMMA
       TempNode{length(TempNode)+1} = Fem.Node;
    end
    
    if ~Fem.SolverStartMMA
        Fem = FiniteElementDataLog(Fem);
    end
    
    if Fem.SolverPlot || (Fem.SolverStartMMA && Fem.SolverPlot)
        Fem.show(Fem.SolverPlotType); 
        drawnow;
    end
    
    if ~Fem.Nonlinear || Fem.Time >= 1
        break; 
    end
end

end
%-------------------------------------------- NL dynamic newmark-beta solve
function Fem = simulate(Fem)
gam               = 1/4; % its not recommended to change these  (1/2)
bet               = 1/2; % its not recommended to change these  (1/4)
Fem.Log           = [];    
Fem.TimeDelta     = Fem.TimeEnd;
Fem.BisectCounter = 1;
Fem.BisectLimit   = 5e3;
Fem.Increment     = 1;
Fem.Iteration     = 1;
Fem.Time          = 0;
Fem.TimeDelta     = 0;
Fem.EndIncrement  = false;
Fem.Convergence   = true;
Fem.TimeStep      = Fem.TimeStep0 + 1e-6;
Fem.Density       = clamp(Fem.Density,Fem.Ersatz,1);
Fem.BisectLimit   = 500;

% if Fem.ResidualNorm > 1e-3
%     Fem.ResidualNorm  = 1e-5; 
%     Fem.StressNorm    = 1e-9; 
%     Fem.DisplaceNorm  = 1e-9; 
% end

if norm(Fem.Utmp) == 0
    Fem.Node      = Fem.Node0;
    Fem.Utmp      = zeros(Fem.Dim*Fem.NNode,1);
else
    Fem.GravityRamp = false;
end

Fem.dUtmp         = zeros(Fem.Dim*Fem.NNode,1);
Fem.ddUtmp        = zeros(Fem.Dim*Fem.NNode,1);
Fem.Ftmp          = zeros(Fem.Dim*Fem.NNode,1);

if Fem.ShowProcess
    showInformation(Fem,'NonlinearFem');
end

if Fem.SolverPlot || (Fem.SolverStartMMA && Fem.SolverPlot)
    figure(101); 
    Fem.show('0');
end

% pre-build
Fem.LoadingFactor = sigmoid(Fem.Time/(5*Fem.TimeStep));
   
% assemble global system matrices
Fem = AssembleGlobalSystem(Fem);

% assemble force/boundary conditions
Fem = AssembleBoundaryConditions(Fem);

% get free DOFs
FreeDofs = GetFreeDofs(Fem);

M  = Fem.MassMatrix(FreeDofs,FreeDofs);
C  = Fem.DampingMatrix(FreeDofs,FreeDofs);
Fi = Fem.fInternal(FreeDofs);
Fe = Fem.fExternal(FreeDofs);
dt = Fem.TimeStep;
  
% precompute acceleration zero
if norm(Fem.Utmp) == 0
    Fem.ddUtmp(FreeDofs) = (M)\(Fe - C*Fem.dUtmp(FreeDofs) - Fi);
else
    % update mesh class
    Fem.Node    = UpdateNode(Fem,Fem.Utmp);
    Fem.Center  = UpdateCenter(Fem);
    Fem.ddUtmp(FreeDofs) = (M)\(Fem.LoadingFactor*Fe ...
        - C*Fem.dUtmp(FreeDofs) - Fem.LoadingFactor*Fi);
end

while true 
    
    [Fem, Terminate] = SolverTimer(Fem);
    
    Fem.Divergence = 1;
    
    if Terminate
        break; 
    end
    
    flag           = 0;
    Singular       = false;
    Fem.Divergence = 0;
    Fem.Utmp_      = Fem.Utmp;
    Fem.dUtmp_     = Fem.dUtmp;
    Fem.ddUtmp_    = Fem.ddUtmp;
   
    % reset iteration
    Fem.Iteration = 1;
    
    % update auxilary flow
    if ~isempty(Fem.Flow)
        Fem.Log.AUX.z = UpdateAuxiliaryFlow(Fem);
    end
    
    while flag == 0
        
        % get free DOFs
        FreeDofs = GetFreeDofs(Fem);
        
        % update intermediate displacement
        ddDelta  = Fem.ddUtmp;

        if Fem.Iteration == 1
            Fem.dUtmp  = Fem.dUtmp_ + ... (1/dt)*(Fem.Utmp - Fem.Utmp_) + ...
                (1-gam)*dt*Fem.ddUtmp_ + gam*dt*ddDelta;
            Fem.Utmp   = Fem.Utmp_ + dt*Fem.dUtmp_ + ...
                (0.5-bet)*dt*dt*Fem.ddUtmp_ + dt*dt*bet*ddDelta;
        end
        
        % compute loading factor (converges to 1 quickly).
        Fem.LoadingFactor = sigmoid(Fem.Time/(5*Fem.TimeStep));
        
        % assemble global system matrices
        Fem = AssembleGlobalSystem(Fem);      
        Fem = AssembleBoundaryConditions(Fem);
        
        % compute tangent stiffness and internal
        M  = Fem.MassMatrix(FreeDofs,FreeDofs);
        R  = Fem.DampingMatrix(FreeDofs,FreeDofs);
        Kt = Fem.TangentStiffness(FreeDofs,FreeDofs);
       
        Fi = Fem.fInternal(FreeDofs);
        Fe = Fem.fExternal(FreeDofs); 
        Fu = Fem.fInput(FreeDofs);
        
        % newmark matrices A        
        A = bet*dt*dt*Kt + gam*dt*R + M;  
        b = Fe + Fu - Fi - R*Fem.dUtmp(FreeDofs);
        
        % solve for acceleration    
        DeltaU = A\b;
        %[DeltaU,~] = gmres(A,b,[],Fem.ResidualNorm,100);
        %[L,D,P] = ldl(A,'vector');
        %DeltaU  = sparse(P,1,(L'\(D\(L\b(P))))); % thanks Ondrej ;)
 
        ddDelta(FreeDofs,1) = ddDelta(FreeDofs,1) + DeltaU(:,1);
        
        Fem.SolverResidual = vappend(Fem.SolverResidual,...
            norm((ddDelta)));
    
        Fem.SolverVonMises = vappend(Fem.SolverVonMises,...
            norm(Fem.s(:,1)));
        
        Fem.SolverDisplace = vappend(Fem.SolverDisplace,...
            norm(b));
        
        % check convergence of solution
        [flag, Fem] = CheckConvergence(Fem,Singular); 
        
        % add iteration
        Fem.Iteration = Fem.Iteration + 1;
        
        % update mesh 
        if flag == 0 && Fem.Nonlinear

            Fem.ddUtmp          = ddDelta;  
            Fem.dUtmp(FreeDofs) = Fem.dUtmp(FreeDofs) + gam*dt*DeltaU;
            Fem.Utmp(FreeDofs)  = Fem.Utmp(FreeDofs)  + bet*dt*dt*DeltaU;
            
        elseif flag == 1 && Fem.Nonlinear
            
            beta = Fem.LoadingFactor;
            
            dq = (1/dt)*(Fem.Utmp - Fem.Utmp_); 
            Fem.Kinetic    = 0.5*beta*dq(FreeDofs).'*M*dq(FreeDofs);
            Fem.PotentialF = Fem.PotentialF + dt*dq(FreeDofs).'*Fu;
            
            Fem.ddUtmp = ddDelta;
            Fem.dUtmp = Fem.dUtmp_ + dt*dt*(1.0-gam)*Fem.ddUtmp_ ...
                + dt*dt*gam*ddDelta;
            
            Fem.Utmp   = Fem.Utmp_ + dt*Fem.dUtmp_ + ...
                (0.5-bet)*dt*dt*Fem.ddUtmp_ + dt*dt*bet*ddDelta;
            
            Fem.Node    = UpdateNode(Fem,Fem.Utmp);
            Fem.Center  = UpdateCenter(Fem);
        end
        
        % update mesh class
        Fem.Mesh = Fem.Mesh.set('Node',Fem.Node);
        Fem.Mesh = Fem.Mesh.set('Center',Fem.Center);
        
    end
    
    if ~Fem.SolverStartMMA
        Fem = AssembleGlobalSystem(Fem);
        Fem = FiniteElementDataLog(Fem);
    end
    
    if Fem.SolverPlot || (Fem.SolverStartMMA && Fem.SolverPlot)
        Fem.show(Fem.SolverPlotType); 
        if ~isempty(Fem.Tendon)
           Fem.show('Fvec'); 
        end
        drawnow;
    end
   
end

end
%---------------------------------------------- topology optimization solve
function Fem = optimize(Fem)
 
if Fem.ShowProcess
    showInformation(Fem,'TopologyOptimization');
end    
    
Fem.SolverPlot     = false;
Fem.IterationMMA   = 0;
Fem.SolverStartMMA = true;
flag               = true;
Visual             = 'ISO';
PenalUpdate        = round(Fem.MaxIterationMMA/5);
Fem.SpatialFilter  = GenerateRadialFilter(Fem,Fem.FilterRadius);

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
    [Fem,ZNew] = UpdateSchemeMMA(Fem,f,dfdz,g,dgdz);
    
    % determine material change
    Fem.Change  = clamp(ZNew - Fem.Density,-Fem.ChangeMax,Fem.ChangeMax);
    Fem.Density = Fem.Density + Fem.Change;
    
    % evaluate fitness
    Fem.Objective  = vappend(Fem.Objective,f);       
    Fem.Constraint = vappend(Fem.Constraint,g);
    
    % check convergence
    [flag,Fem] = CheckConvergenceOpt(Fem);
    
    % draw visual
    if mod(Fem.IterationMMA,10) == 0
        Fem.show(Visual); drawnow;
        colormap(turbo);
    end
    
    if mod(Fem.IterationMMA,PenalUpdate) == 0
       Fem.Penal = clamp(Fem.Penal + 1,1,5);
    end
end

Fem.SolverStartMMA = false;

end
%-------------------------------------------------- compute system matrices
function Fem = compute(Fem,Q,varargin)
    
% get free-dofs    
[dofa,Ia]  = GetFreeDofs(Fem);
Fem.Utmp  = zeros(Fem.Dim*Fem.NNode,1);
Fem.dUtmp = zeros(Fem.Dim*Fem.NNode,1);

if isempty(varargin)
    dQ = Q*0;
end

if size(Q,1) == numel(dofa)
    Fem.Utmp(dofa)  = Q; 
    Fem.dUtmp(dofa) = dQ;
else
    Fem.Utmp(:)  = Q; 
    Fem.dUtmp(:) = dQ;
end

Fem = AssembleGlobalSystem(Fem,1);
Fem = AssembleBoundaryConditions(Fem);

% update nodal;
Fem.Node   = UpdateNode(Fem,Fem.Utmp);
Fem.Center = UpdateCenter(Fem);

Fem.Log.EL.M   = Fem.MassMatrix;
Fem.Log.EL.R   = Fem.DampingMatrix;
Fem.Log.EL.Kt  = Fem.TangentStiffness;
Fem.Log.EL.K   = Fem.Stiffness;
Fem.Log.EL.tau = Fem.fInput;
Fem.Log.EL.fg  = Fem.fExternal;
Fem.Log.EL.fe  = Fem.fInternal;
Fem.Log.dofa   = Ia.';

end
%--------------------------------------------------------------- find nodes
function NodeList = FindNodes(Fem,varargin)
NodeList = FindNode(Fem.Node0,varargin{1:end});
end
%------------------------------------------------------------ find elements
function ElementList = FindElements(Fem,varargin)
ElementList = FindNode(Fem.Mesh.Center,varargin{1:end});
end
%------------------------------------------------------------ find elements
function NodeList = FindEdges(Fem,varargin)
NodeList = FindEdge(Fem.Mesh,varargin{1:end});
end
%---------------------------------------------------------- add constraints
function Fem = AddConstraint(Fem,varargin)   
for ii = 1:3:length(varargin)
  if size(varargin{ii+2},2) == 2 || size(varargin{ii+2},2) == 3 || ...
          isa(varargin{ii+2},'function_handle') || isa(varargin{ii+2},'double')
      if strcmp(varargin{ii},'PressureCell')
         Fem.VolumetricPressure = true; 
         Fem.(varargin{ii}) = [varargin{ii+1},repmat(varargin{ii+2},...
         [length(varargin{ii+1}),1])];
      elseif strcmp(varargin{ii},'Displace')
          Fem.PrescribedDisplacement = true;
          BC = [varargin{ii+1},repmat(transpose(varargin{ii+2}(:)),...
          [length(varargin{ii+1}),1])];
          Fem.Load = [Fem.Load;BC];
      elseif strcmp(varargin{ii},'Contact')
          Fem.(varargin{ii}) = {varargin{ii+1},varargin{ii+2}};
      elseif strcmp(varargin{ii},'Pressure')
          if ~isa(varargin{ii+2},'function_handle')
            Fem.(varargin{ii}) = [Fem.(varargin{ii}) ; ...
                varargin(ii+1),repmat(varargin{ii+2},[length(varargin{ii+1}),1])];
          else
            f = varargin{ii+2};
            Fem.(varargin{ii}) = [Fem.(varargin{ii}); {varargin{ii+1}},{f}];  
          end
      elseif strcmp(varargin{ii},'Gravity')
          Fem.(varargin{ii}) = [varargin{ii+2}(:)];
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
%---------------------------------- export triangular mesh from iso-surface
function msh = exportMesh(Fem,varargin)
    
    ISO = varargin{1};
    [~,I,UxUy] = showISO(Fem,ISO);
      
    B = UxUy;
    Xscale = (B(2)-B(1))/size(I.CData,2);
    Yscale = (B(4)-B(3))/size(I.CData,1);
    
    simplify_tol = varargin{2};
    
    img = I.CData >= 25;
    img = fliplr(img.');
    
    bnd = bwboundaries(img);
    
    c_cell0 = {};
    c_cell = {};
    
    for ii=1:length(bnd)
        bnd_tmp = bnd{ii};
        assert(all(bnd_tmp(1,:)==bnd_tmp(end,:)),'contour is not closed');
        c_cell0{ii} = bnd_tmp;
    end
    
    for ii=1:length(c_cell0)
        c_tmp = c_cell0{ii};
        c_red = decimatePoly(c_tmp,[simplify_tol, 2],false);
        if (nnz(c_red(:,1))>0)&&(nnz(c_red(:,2))>0)
            c_cell{end+1,1} = [Xscale*c_red(:,1), (Yscale)*c_red(:,2)];
        end
    end
    
    % create the 2d triangulation
    H = varargin{3};
    Tesselation = triangulationCreate(c_cell, H(1), H(2), H(3),'linear');
    
    msh = Mesh(Tesselation.Nodes.',Tesselation.Elements.');
    msh = msh.generate();
   
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
cla; axis equal;     
axis off; hold on;

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
axis equal
axis off;

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
    elseif strcmp(varargin{1},'Sdf')
        Pc = Fem.Mesh.get('Center');
        sdf = varargin{2};
        D = sdf(Pc);
        Z = ones(Fem.NElem,1); 
        Z(D(:,end) <= 0) = 0;
        Fem.Density = Z;
    else
        Fem.Density = InitialDesign(Fem,{[1,1],1});
    end
    
    Fem.Density = Fem.SpatialFilter*(1.25*Fem.Density*...
        clamp(Fem.VolumeInfill,0.2,1));
end
%----------------------------------------------------- form smooth topology
function Fem = former(Fem)
    
if Fem.SolverStartMMA && ~Fem.MovieStart
    Res = 100;
    Thickness = 0.2*(Fem.BdBox(2)-Fem.BdBox(1));
else
    Res = 300;
    Thickness = 0.2*(Fem.BdBox(2)-Fem.BdBox(1));
end
    
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

if isempty(Fem.Crop)
    B = Fem.BdBox; 
else
    B = Fem.Crop;
end

x = verts(:,1); 
y = verts(:,2);
dx = 0.01*(B(2) - B(1));
dy = 0.01*(B(4) - B(3));
xq = linspace(B(1)+dx,B(2)-dx,Res); 
yq = linspace(B(3)+dy,B(4)-dy,Res);

[xxq,yyq] = meshgrid(xq,yq);
P = griddata(x,y,V,xxq,yyq);
vrt = [xxq(:),yyq(:)];
Dist = Fem.Mesh.SDF(vrt); 
Dist = Dist(:,end);

tol = 0.1*sqrt((B(2)-B(1))*(B(4)-B(3))/length(vrt));

P = P(:); 
P(Dist > tol) = 0; 
P = (reshape(P,Res,Res)); 
V = cat(3,xxq,yyq,P);

if Fem.VolumetricPressure
    P = griddata(x,y,V2,xxq,yyq);
    Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); 
    Dist = Dist(:,end);
    
    P = P(:); 
    P(Dist > tol) = 0; 
    P(isnan(P))=0;
    P = (reshape(P,Res,Res));
    
    GapFill = cat(3,xxq,yyq,P);
    
    P = griddata(x,y,E,xxq,yyq);
    Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); 
    Dist = Dist(:,end);
    
    P = P(:);  
    P(isnan(P))=0;
    P = (reshape(P,Res,Res));
    
    EdgeFill = cat(3,xxq,yyq,P);
end

SDF0 = V(:,:,3);
SDF0 = GaussianFilter(SDF0,5);
SDF  = repmat(SDF0,[1 1 Layers]);

if Fem.VolumetricPressure
    V       = GapFill;
    SDF2    = V(:,:,3);
    SDF2    = GaussianFilter(SDF2,10);
    ZFiller = repmat(SDF2,[1 1 Patch]);
    
    for ii = 1:Patch
        ZFiller(:,:,ii) = lerp(SDF2+0.01*ii,SDF0,...
                               cos((1/Patch)*pi));
    end
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
function [Fem,I,UxUy] = showISO(Fem,varargin)
V = Fem.Topology;
X = Fem.TopologyGridX; 
Y = Fem.TopologyGridY; 
a = size(V);

if nargin < 3
    depth = floor(a(3)/2);
else
    depth = a(3)*varargin{2}; 
end

X    = X(:,:,depth);
Y    = Y(:,:,depth);
SDF0 = V(:,:,depth);

SDF = GenerateCell(Fem,SDF0);

if ~isempty(Fem.CellRepetion)
    Rpt = Fem.CellRepetion;
    SDF = repmat(SDF,Rpt(2),Rpt(1));
else
    Rpt = [1, 1];
end

if ~isempty(Fem.ReflectionPlane)
    Rp = Fem.ReflectionPlane;
    
    if Rp(1) == 1
        Xf = flip(X) + (max(max(max(X))) - min(min(min(X))));
        X  = vertcat(Xf,X);
    end
    
    if Rp(1) == -1
        Xf = flip(X) - (max(max(max(X))) - min(min(min(X))));
        X  = horzcat(X,Xf);
    end
    
    if Rp(2) == -1
        Yf = flip(Y) - (max(max(max(Y))) - min(min(min(Y))));
        Y  = horzcat(Y,Yf);
    end
    
    if Rp(2) == 1
        Yf = flip(Y) + (max(max(max(Y))) - min(min(min(Y))));
        Y  = horzcat(Yf,Y);
    end
end

scaleX = 1; scaleY = 1;
if ~isempty(Fem.Repeat)
    
    SDF0     = SDF;
    Instruct = Fem.Repeat;
    
    for ii = 1:length(Instruct)
        
        if Instruct(ii) == 1
            scaleX = scaleX + 1;
            SDF    = cat(2,SDF,SDF0);
        end
        
        if Instruct(ii) == 2
            scaleY = scaleY*2;
            SDF    = cat(1,SDF,SDF);
        end
    end
end

cla;
Uxx  = scaleX*Rpt(1)*[min(min(min(X))) max(max(max(X)))];
Uyy  = scaleY*Rpt(2)*[min(min(min(Y))) max(max(max(Y)))];
UxUy = [0, Uxx(2)-Uxx(1), 0, Uyy(2)-Uyy(1)];

I = GaussianFilter((SDF >= varargin{1})*255,3);

I = image(rescale(Uxx),...
    ((max(Uyy)-min(Uyy))/(max(Uxx) - min(Uxx)))*rescale(Uyy),I);

axis equal; axis off; 
colormap(Fem.ColormapOpt); 
caxis([0 1]);
background('w');
end
end

methods (Access = private)
%%///////////////////////////////////////////////////////// FINITE ELEMENTS
%------------------------------------------------ convert mesh to fem class
function Fem = SetupFiniteElement(Fem)

Fem.ElemNDof      = Fem.Dim*cellfun(@length,Fem.Element);
Fem.SpatialFilter = GenerateRadialFilter(Fem,Fem.FilterRadius);

if (~Fem.AssembledSystem && ~Fem.Nonlinear) || Fem.Nonlinear
Fem.i = zeros(sum(Fem.ElemNDof.^2),1);
Fem.j  = Fem.i; Fem.e  = Fem.i; Fem.fi = Fem.i; Fem.k  = Fem.i; 
Fem.m  = Fem.i; Fem.c  = Fem.i; Fem.t  = Fem.i; Fem.fb = Fem.i; 
Fem.ft = Fem.i; Fem.v  = Fem.l;
Fem.s  = zeros(Fem.NNode,6); 
Fem.s  = zeros(Fem.NNode,3);
Fem.l  = zeros(Fem.NNode,1); 

Fem.Potential  = 0;
Fem.PotentialG = 0;
Fem.PotentialF = 0;
Fem.Kinetic    = 0;
end

switch Fem.Solver
    case('Basic');    Fem.SolverId = 1; % mrdivide
    case('mr');       Fem.SolverId = 1; % mrdivide
    case('lu');       Fem.SolverId = 2; % LU-decomp min eig-value   
    case('gmres');    Fem.SolverId = 3; % Generalized Minimum Residual           
    otherwise;        Fem.SolverId = 1;
end

end
%------------------------------------ assemble global finite-element system 
function Fem = AssembleGlobalSystem(Fem,ForceBuild)
    
if nargin < 2
    ForceBuild = false; 
end

% evaluate shape-functions at nodal locations
if (Fem.Iteration == 1 && Fem.Increment == 1)
    tab = struct; tab.Element = Fem.Element;
    
    switch(Fem.Mesh.Type)
        case('C2PX'), tab = TabulateShapeFunctions(tab);
        case('C2T3'), tab = TabulateShapeFunctions(tab);
        case('C2Q4'), tab = TabulateShapeFunctions(tab);
        case('C3H8'), tab = TabulateShapeFunctionsC3H8(tab);
        case('C3T4'), tab = TabulateShapeFunctionsC3T4(tab);
        otherwise,    tab = TabulateShapeFunctions(tab);
    end
    
    Fem.ShapeFnc = tab.ShapeFnc;
end

if Fem.VolumetricPressure
    beta = Fem.OptFactor*Fem.LoadingFactor;
    dV   = beta*Fem.PressureCell(:,2)*ones(Fem.NElem,1);
elseif ~isempty(Fem.Contraction)
    beta = Fem.OptFactor*Fem.LoadingFactor;
    
    if isa(Fem.Contraction,'funtion_handle')
       dV = Fem.Contraction(Fem);
    else
       dV = beta*Fem.Contraction;
    end
else 
    dV   = 0*ones(Fem.NElem,1);
end

beta = Fem.LoadingFactor;
[E,~,V] = MaterialField(Fem);

Fem.Potential  = 0;
Fem.PotentialG = 0;

index    = 0; 
subindex = 0;

if strcmp(Fem.Material.Type,'NeoHookean')
    LocalType = 1;
elseif strcmp(Fem.Material.Type,'Yeoh')
    LocalType = 2;
elseif strcmp(Fem.Material.Type,'Mooney')
    LocalType = 3;
end

if (~Fem.AssembledSystem && ~Fem.Nonlinear) ...
        || Fem.Nonlinear || ForceBuild
    
for el = 1:Fem.NElem
   
    NDof = Fem.ElemNDof(el);
    if Fem.Dim == 2
        eDof = reshape([2*Fem.Element{el}-1;
            2*Fem.Element{el}],NDof,1);
    else
        eDof = reshape([3*Fem.Element{el}-2;
            3*Fem.Element{el}-1;
            3*Fem.Element{el}],NDof,1);
    end

    if LocalType == 1
        nn = length(Fem.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsNHFast(Fem.Element{el},eDof,dV(el),V(el),...
            Fem.Dim,Fem.Node0,...
            Fem.ShapeFnc{nn}.N,...
            Fem.ShapeFnc{nn}.dNdxi,...
            Fem.ShapeFnc{nn}.W,...
            Fem.Utmp, Fem.dUtmp,...
            Fem.Material.Rho, Fem.Material.Zeta,...
            Fem.Gravity,...
            Fem.Material.Mu,...
            Fem.Material.Lambda);
    elseif LocalType == 2
        nn = length(Fem.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsYHFast(Fem.Element{el},eDof,dV(el),V(el),...
            Fem.Dim,Fem.Node0,...
            Fem.ShapeFnc{nn}.N,...
            Fem.ShapeFnc{nn}.dNdxi,...
            Fem.ShapeFnc{nn}.W,...
            Fem.Utmp, Fem.dUtmp,...
            Fem.Material.Rho, Fem.Material.Zeta,...
            Fem.Gravity,...
            [Fem.Material.C1,Fem.Material.C2,Fem.Material.C3],...
            [Fem.Material.D1,Fem.Material.D2,Fem.Material.D3]);    
    elseif LocalType == 3
        nn = length(Fem.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsMNFast(Fem.Element{el},eDof,dV(el),V(el),...
            Fem.Dim,Fem.Node0,...
            Fem.ShapeFnc{nn}.N,...
            Fem.ShapeFnc{nn}.dNdxi,...
            Fem.ShapeFnc{nn}.W,...
            Fem.Utmp, Fem.dUtmp,...
            Fem.Material.Rho, Fem.Material.Zeta,...
            Fem.Gravity,...
            [Fem.Material.C10,Fem.Material.C01],...
            Fem.Material.K);       
%     else
%         [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
%             Locals(Fem,Fem.Element{el},eDof,dV,V(el));
    end
     
    ind1 = index+1:index+NDof^2;               % matrix indexing
    ind2 = index+1:index+NDof;                 % vector indexing
    ind3 = subindex+1:subindex+NDof/Fem.Dim;   % dimension indexing

    I = repmat(eDof,1,NDof); J = I';
    Fem.e(ind1)  = el;
    Fem.i(ind1)  = I(:);
    Fem.j(ind1)  = J(:);
    Fem.m(ind1)  = Me(:);
    Fem.c(ind1)  = Ce(:);
    Fem.k(ind1)  = Ke(:);
    Fem.t(ind1)  = Kte(:);
    Fem.fi(ind2) = Fe(:);
    Fem.fb(ind2) = Fb(:);
    Fem.ft(ind2) = Te(:);    
    
    Fem.s(ind3,1) = E(el)*Svme(:);
    Fem.s(ind3,2) = E(el)*SS(:,1);
    Fem.s(ind3,3) = E(el)*SS(:,2);
    Fem.s(ind3,4) = E(el)*SS(:,4);
    Fem.p(ind3,1) = E(el)*EE(:,1);
    Fem.p(ind3,2) = E(el)*EE(:,2);
    Fem.p(ind3,3) = E(el)*EE(:,4);
    Fem.l(ind3)   = Fem.Element{el}(:);
    
    Fem.Potential  = Fem.Potential  + beta*Ve;
    Fem.PotentialG = Fem.PotentialG + beta*Vge;
    Fem.v(ind3)    = Qe(:);
    
    if Fem.Iteration == 1
        Fem.ElemRot{el} = Re;
        Fem.ElemStr{el} = Ue;
    end
    
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

% build global matrix of additional springs
spMat = spdiags(sp(:),0,Fem.Dim*Fem.NNode,Fem.Dim*Fem.NNode);

% build output vector of nodal displacements
NOutput = size(Fem.Output,1);
L       = sparse(Fem.Dim*Fem.NNode,1);

if NOutput ~=0
    L(2*Fem.Output(1:NOutput,1)-1) = Fem.Output(1:NOutput,2);
    L(2*Fem.Output(1:NOutput,1))   = Fem.Output(1:NOutput,3);
end
    
K   = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);  % init global stiffness
Ktr = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);  % init global tangent stiffness
Fe  = sparse(Fem.i,1,E(Fem.e).*Fem.fi);     % init elastic force vector
Fb  = sparse(Fem.i,1,E(Fem.e).*Fem.fb);     % init body force vector
Ft  = sparse(Fem.NNode*Fem.Dim,1);          % init contraction force
F   = sparse(Fem.Dim*Fem.NNode,1);          % init global force vector
M   = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.m);  % init mass matrix
C   = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.c);  % init dampings matrix

if Fem.Nonlinear
    beta = Fem.OptFactor*Fem.LoadingFactor;
else
    beta = 1;
end

if ~isempty(Fem.Load) && ~Fem.PrescribedDisplacement
    NLoad = size(Fem.Load,1);
    
    for ii = 1:Fem.Dim
        F(Fem.Dim*Fem.Load(1:NLoad,1)+(ii-Fem.Dim),1) = ...
            beta*Fem.Load(1:NLoad,1+ii);
    end
end

if ~isempty(Fem.Tendon) && ~Fem.PrescribedDisplacement
    F_     = F*0;
    List     = 1:Fem.NElem;
    N2F      = Fem.Mesh.NodeToFace;
    NTendons = size(Fem.Tendon,1);
    
    for jj = 1:NTendons
        
        el  = List(N2F(Fem.Tendon(jj,1),:) > 0);
        
        % compute mean rotation matrix between elements
        Rot = zeros(3);
        for kk = 1:numel(el)
            Rot = Rot + (1/numel(el))*Fem.ElemRot{el(kk)};
        end
        
        % ensure orthogonal
        [Ur,~,Vr] = svd(Rot);
        Rot = (Ur*Vr.');
        
        if Fem.Dim == 2
            Ft = Rot(1:2,1:2)*Fem.Tendon(jj,2:end).';
        else
            Ft = Rot*Fem.Tendon(jj,2:end).';
        end
        
        for ii = 1:Fem.Dim
            F_(Fem.Dim*Fem.Tendon(jj,1)+(ii-Fem.Dim),1) = ...
                beta*Ft(ii);
        end
    end
    
    F = F + F_;
end

if ~isempty(Fem.Contact) && Fem.PrescribedDisplacement == true 
    
    Y    = Fem.Node;
    SDF  = Fem.Contact{1};
    Move = Fem.Contact{2};
    
    Y(:,1) = Y(:,1) - beta*Move(1);
    Y(:,2) = Y(:,2) - beta*Move(2);
    
    d = SDF(Y); I = find((d(:,end))<0);   
    eps = 1e-5;
    n1 = (SDF(Y(I,:)+repmat([eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
    n2 = (SDF(Y(I,:)+repmat([0,eps],size(Y(I,:),1),1))-d(I,end))/eps;

    F_     = F*0;
    F_(2*I(1:size(Y(I,:),1),1)-1,1) = -1.99*d(I,end).*n1(:,end); 
    F_(2*I(1:size(Y(I,:),1),1),1)   = -1.99*d(I,end).*n2(:,end);
    
    pDof = zeros(2*length(I),1);
    
    for kk = 1:length(I)
        pDof(2*kk-1) = 2*I(kk)-1; 
        pDof(2*kk)   = 2*I(kk);
    end

    F = F + F_;
end

if ~isempty(Fem.Contact) && Fem.PrescribedDisplacement == false 
    
    F_ = F*0;
    SDF  = Fem.Contact{1};
    Move = Fem.Contact{2};
    Emod = Fem.Material.getModulus();
    Cmod = Fem.Material.getContactFriction();
    
    Y0 = Fem.Node0;
    
    [Ux,Uy]   = DisplacementField(Fem,Fem.Utmp);
    [dUx,dUy] = DisplacementField(Fem,Fem.dUtmp);
    
    V  = sqrt(dUx.^2 + dUy.^2);
    Vc = mean(V);
    
    Y(:,1) = Y0(:,1) + Ux - beta*Move(1);
    Y(:,2) = Y0(:,2) + Uy - beta*Move(2);
    
    if Fem.Dim == 3
        Y(:,3) = Y0(:,3) + Uy - beta*Move(3);
    end
    
    eps = 1e-2;
    
    d = SDF(Y); 
    I = find((d(:,end))<eps);
    
    if ~isempty(I)
        if Fem.Dim == 3
            n1 = (SDF(Y(I,:)+repmat([eps,0,0],size(Y(I,:),1),1))-d(I,end))/eps;
            n2 = (SDF(Y(I,:)+repmat([0,eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
            n3 = (SDF(Y(I,:)+repmat([0,0,eps],size(Y(I,:),1),1))-d(I,end))/eps;
            
            Ux = (d(I,end)).*n1(:,end);
            Uy = (d(I,end)).*n2(:,end);
            Uz = (d(I,end)).*n3(:,end);
            vUx = (dUx(I,end)).*n1(:,end);
            vUy = (dUy(I,end)).*n2(:,end);
            vUz = (dUz(I,end)).*n3(:,end);
            
            N(:,1) = n1(:,end);
            N(:,2) = n2(:,end);
            N(:,3) = n3(:,end);
            
            N = N./sqrt((sum((N.^2),2)));
            T = ([0,1,0;-1,0,0;0,0,1]*N.').';
            
            Ffric = -0*N.*dot([vUx,vUx*0,vUx*0].',T.').';
            
            F_(3*I(1:size(Y(I,:),1),1)-2,1) = -0.1*Emod*(Ux).^3 - 0*Emod*Ux.*abs(vUx).^2 ...
                + Emod*abs(d(I,end)).*Ffric(:,1);
            
            F_(3*I(1:size(Y(I,:),1),1)-1,1)   = -0.1*Emod*(Uy).^3 - 0*Emod*Uy.*abs(vUy).^2 ...
                + Emod*abs(d(I,end)).*Ffric(:,2);
            
            F_(3*I(1:size(Y(I,:),1),1),1)   = -0.1*Emod*(Uz).^3 - 0*Emod*Uz.*abs(vUz).^2 ...
                + Emod*abs(d(I,end)).*Ffric(:,3);
            
            C(3*I(1:size(Y(I,:),1),1)-1,2*I(1:size(Y(I,:),1),1)-1) = ...
                C(2*I(1:size(Y(I,:),1),1)-1,2*I(1:size(Y(I,:),1),1)-1) + Cmod;
            
            C(3*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) = ...
                C(2*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) + Cmod;
            
            C(3*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) = ...
                C(2*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) + Cmod;
            
        else
            n1 = (SDF(Y(I,:)+repmat([eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
            n2 = (SDF(Y(I,:)+repmat([0,eps],size(Y(I,:),1),1))-d(I,end))/eps;
            
            Ux = (d(I,end)).*n1(:,end);
            Uy = (d(I,end)).*n2(:,end);
            vUx = (dUx(I,end)).*n1(:,end);
            vUy = (dUy(I,end)).*n2(:,end);
            
            N(:,2) = n1(:,end);
            N(:,1) = -n2(:,end);
            T = N./sqrt((sum((N.^2),2)));
            
            Ffric = -0*N.*dot([vUx,vUy].',T.').';
            
            id = I(1:size(Y(I,:),1),1);
            
            % modify forcing
            F_(2*id-1,1) = -0.1*Emod*(Ux).^3 - 0.0*Emod*Ux.*abs(vUx).^2 ...
                + 1e-4*Fem.Material.Zeta*Emod*abs(d(I,end)).*Ffric(:,1);
            F_(2*id,1)   = -0.1*Emod*(Uy).^3 - 0.0*Emod*Uy.*abs(vUy).^2 ...
                + 1e-4*Fem.Material.Zeta*Emod*abs(d(I,end)).*Ffric(:,2);
            
            % modify damping
            C(2*id-1,2*id-1) = C(2*id-1,2*id-1) + Cmod*abs(d(I,end)).^2;
            C(2*id,2*id)     = C(2*id,2*id) + Cmod*abs(d(I,end)).^2;
            
        end
        
        F = F + F_;
    end
end

if Fem.PrescribedDisplacement
    NLoad = size(Fem.Load,1);
    
    if ~isempty(Fem.Load)
        pDof = [];
        
        if (numel(Fem.Load(1,2:end)) < 6 && Fem.Dim == 3) || ...
                (numel(Fem.Load(1,2:end)) < 4 && Fem.Dim == 2)
            for ii = 1:Fem.Dim
                pd = reshape(Fem.Dim*Fem.Load(:,1) + ...
                             (ii-Fem.Dim),NLoad,1);
                F(pd) = Fem.Load(1:NLoad,1+ii);
                if norm(Fem.Load(1:NLoad,1+ii)) > 0
                    pDof = [pDof;pd];
                end
            end
        else
            R  = reshape(Fem.Load(1,2:end),[Fem.Dim,Fem.Dim]);
            N0 = mean(Fem.Node(Fem.Load(:,1),:),1);
            V  = (R*(Fem.Node(Fem.Load(:,1),:)-N0).').';
            
            for ii = 1:Fem.Dim
                pd = reshape(Fem.Dim*Fem.Load(:,1) + ...
                             (ii-Fem.Dim),NLoad,1);
                F(pd) = V(:,ii);
                pDof  = [pDof; pd];
            end
        end
    end
    
    EMod  = Fem.Material.getModulus();
    I0    = eye(Fem.Dim*Fem.NNode,Fem.Dim*Fem.NNode);
   
    if Fem.Nonlinear
        KTtmp = full(Ktr);
        KTtmp(pDof,:)    = 0*I0(pDof,:);
        KTtmp(pDof,pDof) = -EMod*eye(length(pDof));
        Ktr = KTtmp;

        if isempty(Fem.Contact) 
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
        
        Ktmp         = full(K); 
        Ktmp(pDof,:) = I0(pDof,:);
        Ktmp(:,pDof) = I0(:,pDof);
        
        F = -beta*K*F;
        K = Ktmp;
    end
end

if ~isempty(Fem.Contraction)
    %Fcntr = sparse(Fem.i,Fem.e,Fem.ft);
    %F(:,1) = beta*Fcntr*Fem.Contraction(:);
    W  = ones(Fem.NElem,1);
    Ft = beta*sparse(Fem.i,1,W(Fem.e).*Fem.ft);
end

if ~isempty(Fem.Pressure)
    
    Nds      = Fem.Node0;       
    Nds(:,1) = Nds(:,1) + Fem.Utmp(1:2:end-1,1);
    Nds(:,2) = Nds(:,2) + Fem.Utmp(2:2:end,1);
    
    for kk = 1:size(Fem.Pressure,1)
    EdgeList = Fem.Pressure{kk,1};
    
    if ~isa(Fem.Pressure{kk,2},'function_handle')
        Pload = mean(beta*Fem.Pressure{kk,2});
    else     
        Pload = Fem.Pressure{kk,2}(Fem);
    end
    
    for ii = 1:length(EdgeList)
        NodeID = EdgeList{ii};
        V      = Nds(NodeID,:);

        dsx = diff(V(:,1));
        dsy = diff(V(:,2));
        dl  = sqrt(dsx.^2+dsy.^2);
        Nx  = -dsy./dl;
        Ny  = dsx./dl;
        
        S = [NodeID(1:end-1), NodeID(2:end)].';

        for jj = 1:length(NodeID)-1
            F(Fem.Dim*S(1,jj)-1,1) = F(Fem.Dim*S(1,jj)-1,1) ...
                + 0.5*Pload*Nx(jj)*dl(jj);
            F(Fem.Dim*S(2,jj)-1,1) = F(Fem.Dim*S(2,jj)-1,1) ...
                + 0.5*Pload*Nx(jj)*dl(jj);
            F(Fem.Dim*S(1,jj),1) = F(Fem.Dim*S(1,jj),1) ...
                + 0.5*Pload*Ny(jj)*dl(jj);
            F(Fem.Dim*S(2,jj),1) = F(Fem.Dim*S(2,jj),1) ...
                + 0.5*Pload*Ny(jj)*dl(jj);
           
        end
        
    end
    
    end
end

if Fem.VolumetricPressure
    
    Pc = Fem.Mesh.get('Center');
    W  = ones(Fem.NElem,1);
    
    id     = FindElements(Fem,'FloodFill',Fem,Fem.Density);
    idnull = setdiff(1:Fem.NElem,id);
    
    W(idnull) = 0;
    
    Ft = beta*sparse(Fem.i,1,W(Fem.e).*Fem.ft);
    
    if (Fem.SolverStartMMA || Fem.Nonlinear) 
        z0 = Fem.Density;
        dz = Fem.ChangeMax;
        Fem.dFdE = zeros(2*Fem.NNode,Fem.NElem);
        bnd = boundary(Pc(id,1),Pc(id,2));
%         idbound = id(bnd);
%         for ii = idbound
%             ze = z0;
%             ze(ii) = ze(ii) + dz;
%             id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
%             idnull = setdiff(1:Fem.NElem,id);
%             W = ones(Fem.NElem,1);
%             W(idnull) = 0;
%             SS = dot(W(Fem.e),Fem.ft);
%             fe = sparse(Fem.i,1,SS);
%             Fem.dFdE(:,ii) = ((fe - Ft)/(dz));
%         end
%         Fem.Density = z0;
    end
end


if Fem.GravityRamp
    Fem.fExternal    = beta*Fb;
else
    Fem.fExternal    = Fb;
end

Fem.fInternal        = Fe; 
Fem.fInput           = F + Ft;
Fem.Stiffness        = K + spMat;
Fem.TangentStiffness = Ktr + spMat;
Fem.Residual         = Fem.fInternal - Fem.fExternal - Fem.fInput;
Fem.OutputVector     = L;
Fem.MassMatrix       = M;
Fem.DampingMatrix    = C;

end
%---------------------------------------------- get free degrees-of-freedom
function [FreeDofs, Ia] = GetFreeDofs(Fem)
FixedDofs = [];
NSupp     = size(Fem.Support,1);

if Fem.Dim == 2 && NSupp > 0
    FixedDofs = [Fem.Support(1:NSupp,2).*(2*Fem.Support(1:NSupp,1)-1);
                 Fem.Support(1:NSupp,3).*(2*Fem.Support(1:NSupp,1))];
elseif Fem.Dim == 3 && NSupp > 0
    FixedDofs = [Fem.Support(1:NSupp,2).*(3*Fem.Support(1:NSupp,1)-2);
                 Fem.Support(1:NSupp,3).*(3*Fem.Support(1:NSupp,1)-1);
                 Fem.Support(1:NSupp,4).*(3*Fem.Support(1:NSupp,1))];
end

FixedDofs  = unique(FixedDofs(FixedDofs>0));
AllDofs    = 1:Fem.Dim*Fem.NNode;
[FreeDofs] = setdiff(AllDofs,FixedDofs);
[Ia,~] = ismember(AllDofs(:),FreeDofs(:));
end
%--------------------------------------------------- local element matrices
function [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,Tke,Re,Ue] = ...
        Locals(Fem,eNode,eDof,dV,Rb)
% get order
nn = length(eNode);
mm = Fem.Dim;

% get gauss points and weights
W   = Fem.ShapeFnc{nn}.W;
Fe  = zeros(mm*nn,1);
Fb  = zeros(mm*nn,1);
Te  = zeros(mm*nn,1);
Me  = zeros(mm*nn,mm*nn);
Ce  = zeros(mm*nn,mm*nn);
Ke  = zeros(mm*nn,mm*nn);
Kte = zeros(mm*nn,mm*nn);
SGP = zeros(length(W),6);
EGP = zeros(length(W),6);
RRe = zeros(3);
UUe = zeros(3);
Qe  = ones(nn,1);
Ve  = 0;
Vge = 0;

Nshp = length(Fem.ShapeFnc{nn}.N(:,:,1));
NNe  = zeros(length(W),Nshp);

if Fem.Dim == 2
    Et = [dV; dV; 0];
else
    Et = [dV; dV; dV; 0; 0; 0];
end

% get displacement field
Delta  = Fem.Utmp(eDof,:);
dDelta = Fem.dUtmp(eDof,:);

% computation of Piolla stress
PiollaStress = @(x) Fem.Material.PiollaStress(x);

% quadrature loop
for q = 1:length(W)
    
    % extract shape-functions
    N     = Fem.ShapeFnc{nn}.N(:,:,q);
    dNdxi = Fem.ShapeFnc{nn}.dNdxi(:,:,q);
    J0    = Fem.Node0(eNode,:).'*dNdxi;
    dNdx  = dNdxi/J0;
    dJ    = det(J0);
    
    % deformation gradient   
    F = DeformationGradient(Fem,Delta,dNdx);
    
    % polar decompostion
    %[R, Q, ~] = PolarDecomposition(Fem,F);
    [R, Q, ~] = PolarDecompositionFast_mex(F);

    % increase robustness low density in topo
    Q = Rb*(Q-eye(3)) + eye(3);
    
    % reconstruct deformation gradient
    F = R*Q;
    
%     % right cauchy-green strain
%     C = F.'*F;
%     
%     if strcmp(Fem.Type,'PlaneStress') && Fem.Dim < 3
%         C(3,3) = det(F)/(C(1,1)*C(2,2) - C(1,2)*C(2,1));
%     end

    % get internal stress matrix
    [S0, D0, Psi] = PiollaStress(F);
    %[S0, D0, Psi] = PiollaStressNHFast_mex(F,C1,K);
    
    % voigt-notation vectorize
    S = VoightNotation(S0);
    
    % reduced isotropic matrices
    [Se, De, Ge] = IsotropicReduction(Fem,D0,S);
    
    % nonlinear strain-displacement operator
    [Bnl,Bg,NN,tau] = NonlinearStrainOperatorFast_mex(N,dNdx,F);
    
    % local elemental rotation
    RRe = RRe + R/nn;
    UUe = UUe + Q/nn;
    
    % internal force vector
    Fe = Fe + tau*W(q)*Bnl.'*Se*dJ;
    
    % (graviational) body force vector
    Fb = Fb + tau*W(q)*Fem.Material.Rho*(NN.')*Fem.Gravity(:)*dJ;
    
    % lineararized stiffness matrix
    Ke = Ke + tau*W(q)*(Bnl.'*De*Bnl)*dJ;
    
    % tangent stiffness matrix
    Kte = Kte + tau*W(q)*(Bnl.'*De*Bnl + Bg.'*Ge*Bg)*dJ;
    
    % ensure its symmetric
    Kte = 0.5*(Kte + Kte.');
    
    % mass matrix
    Me = Me + tau*W(q)*Fem.Material.Rho*(NN.')*NN*dJ;
    
    % dampings matrix
    Ce = Ce + Fem.Material.Zeta*Me;
    
    % thermal expansion force
    Te = Te + tau*W(q)*Bnl.'*De*Et*dJ;
    
    % elemental potential energy
    Ve = Ve  + tau*W(q)*Psi*dJ;
    
    % graviational energy
    Vge = Vge - tau*W(q)*((N.'*Fem.Node0(eNode,:)) + ...
        (NN*Delta).')*Fem.Material.Rho*Fem.Gravity(:)*dJ;
    
    % lagrangian strain
    Elagran = (1/2)*(F.'*F - eye(3));
    
    % true stress and lagrangian strain
    SGP(q,:) = VoightNotation((1/det(F))*F*S0*(F.'));
    EGP(q,:) = VoightNotation(Elagran);
    
    % construct shape functions
    NNe(((q-1)*Nshp + 1):(q*Nshp)) = N(:).';
end

% compute elemental kinetic energy
Tke = 0.5*(dDelta).'*Me*(dDelta);    

%compute elemental rotation matrix
[Ur,~,Vr] = svd(RRe);
Re = (Ur*Vr.');
Ue = UUe;

SS = NNe.'*SGP;
EE = NNe.'*EGP;
[Svm, ~] = VonMises(SS(:,1),SS(:,2),SS(:,3),...
                    SS(:,4),SS(:,5),SS(:,6));
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
C = F'*F;
[Q0, lambdasquare] = eig(C);

lambda = sqrt(diag((lambdasquare))); 
Uinv   = repmat(1./lambda',size(F,1),1).*Q0*Q0';

R = F*Uinv;
S = R'*F;
V = F*R';
end
%---------------------------------------------------------- strain operator
function [Bn,Bg,NN,tau] = NonlinearStrainOperator(~,N,dNdx,F)
nn = length(N);
mm = size(dNdx,2);
zz = mm*nn;

NN = zeros(mm,zz);

id1 = 1:mm:zz;
id2 = 2:mm:zz;

NN(1,id1) = N.';
NN(2,id2) = N.';
dNdxX = dNdx(:,1).';
dNdxY = dNdx(:,2).';

if mm == 3
    NN(3,3:mm:zz) = N.';
    dNdxZ = dNdx(:,3).';
end

Bn = zeros((mm-1)*3,zz);
Bg = zeros((mm-1)*4+(mm-2),zz);

if mm == 2 % 2-dimensional
    Bn(1,id1) = dNdxX*F(1,1);
    Bn(1,id2) = dNdxX*F(2,1);
    Bn(2,id1) = dNdxY*F(1,2);
    Bn(2,id2) = dNdxY*F(2,2);
    Bn(3,id1) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(3,id2) = dNdxX*F(2,2) + dNdxY*F(2,1);
    
    Bg(1,id1) = dNdxX;
    Bg(2,id1) = dNdxY;
    Bg(3,id2) = dNdxX;
    Bg(4,id2) = dNdxY;
    
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
if Fem.Dim == 2
    G = zeros(4,4);
    D   = [D0(1,1), D0(1,2),       0;
           D0(2,1), D0(2,2),       0;
                 0,       0, D0(4,4)];
    SIG = [S0(1), S0(4); S0(4), S0(2)];
    S   = [S0(1); S0(2); S0(4)]; 
    
    G(1:2,1:2) = SIG;
    G(3:4,3:4) = SIG;
else
    G = zeros(9,9);
    D = D0;
    S = S0;
    SIG = [S0(1), S0(4), S0(6);
           S0(4), S0(2), S0(5);
           S0(6), S0(5),S0(3)];
       
    G(1:3,1:3) = SIG;
    G(4:6,4:6) = SIG;
    G(7:9,7:9) = SIG;
end

end
%------------------------------------------------------ isotropic reduction
function A = HessianBuild(Fem,Mi,Kt,Rt)
alpha = 1/2;
dt    = Fem.TimeStep;
n     = size(Mi,1);

F = zeros(2*n,2*n);
F(1:n,n+1:2*n)     = Mi;
F(n+1:2*n,1:n)     = -Kt;
F(n+1:2*n,n+1:2*n) = -Rt*Mi;

A = eye(2*n) - alpha*dt*F;

end
%------------------------------------------------------ isotropic reduction
function [flag, Fem] = CheckConvergence(Fem,SingularKt)

CriteriaResidual = (Fem.SolverResidual(end) > Fem.ResidualNorm);

DiffSVM = abs(Fem.SolverVonMises(end-1)/Fem.SolverVonMises(end) - 1);

CriteriaStress = (DiffSVM > Fem.StressNorm);

Criteria = CriteriaResidual && CriteriaStress; %&& CriteriaDisplace;

if Fem.SolverResidual(end,1) > Fem.SolverResidual(end-1,1)
   Fem.Divergence = Fem.Divergence + 1;
end
 
if (Criteria && (Fem.Iteration <= Fem.MaxIteration) && ~SingularKt && ...
         Fem.Divergence < Fem.BisectLimit && Fem.Iteration <= 1)
    flag = 0;
else
    if (Fem.Iteration > Fem.MaxIteration) || SingularKt %|| ...
            %Fem.SolverResidual(end) > (25*Fem.Material.getModulus)
        flag = 2; 
        Fem.Convergence = false;
        Fem.Utmp = Fem.U;
    
    elseif (Fem.Divergence >= Fem.BisectLimit)
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

if Fem.ShowProcess
    ProcessMonitor(Fem);
end

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
    
    if Fem.BisectCounter == 1 && Fem.SolverId ~= 3
        Fem.SolverId = 3;
        if ~Fem.SolverStartMMA
        fprintf('------------------------------------------------------------|\n');
        fprintf(' Switching solver -> Generalized Min. Residual Method       |\n');
        fprintf('------------------------------------------------------------|\n');
        end
        Fem.BisectCounter = Fem.BisectCounter + 1;
        return;
    end
    
    Fem.Utmp = Fem.U;
    Fem.TimeStep = clamp(Fem.TimeStep/2,Fem.TimeStepMin,1);
    Fem.Time = clamp(Fem.Time - Fem.TimeStep,0,1);
    Fem.BisectCounter = Fem.BisectCounter + 1;
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
%-------------------------------------------------- solver logging function
function Fem = FiniteElementDataLog(Fem)
        
% record the stresses for all nodes
Fem.VonMisesNodal = full(sparse(Fem.l,1,Fem.s(:,1))./sparse(Fem.l,1,Fem.v));
Fem.sxxNodal = full(sparse(Fem.l,1,Fem.s(:,2))./sparse(Fem.l,1,Fem.v));
Fem.syyNodal = full(sparse(Fem.l,1,Fem.s(:,3))./sparse(Fem.l,1,Fem.v));
Fem.sxyNodal = full(sparse(Fem.l,1,Fem.s(:,4))./sparse(Fem.l,1,Fem.v));
Fem.exxNodal = full(sparse(Fem.l,1,Fem.p(:,1))./sparse(Fem.l,1,Fem.v));
Fem.eyyNodal = full(sparse(Fem.l,1,Fem.p(:,2))./sparse(Fem.l,1,Fem.v));
Fem.exyNodal = full(sparse(Fem.l,1,Fem.p(:,3))./sparse(Fem.l,1,Fem.v));

% record the forces for all nodes
force        = full(sparse(Fem.i,1,Fem.fi));
Fem.fxNodal  = force(2*(1:Fem.NNode)-1);
Fem.fyNodal  = force(2*(1:Fem.NNode));
Fem.rotNodal = 0*force((1:Fem.NNode));

list = 1:Fem.NElem;
F2V  = Fem.Mesh.FaceToNode;

% for eacth element 
for ii = 1:Fem.NNode
    Rtmp = 0;
    Qtmp = 0;
    ElemFromNode = list(F2V(:,ii) > 0);
    
    for jj = ElemFromNode
        Qtmp = Qtmp + Fem.ElemStr{jj};
        Rtmp = Rtmp + Fem.ElemRot{jj};
    end
    
    [Ur,~,Vr] = svd(Rtmp);
    R{ii} = (Ur*Vr.');
    Q{ii} = Qtmp/numel(ElemFromNode);
end

[ux,uy,un] = DisplacementField(Fem,Fem.Utmp);
[fx,fy,fn] = DisplacementField(Fem,Fem.fInternal);

% construct log file 
if isempty(Fem.Log)
   Fem.Log     = struct;
   Fem.Log.t   = [0;Fem.Time]; 
   Fem.Log.Psi = [0;Fem.Potential];
   Fem.Log.Vg  = [Fem.PotentialG;Fem.PotentialG];
   Fem.Log.Vf  = [Fem.PotentialF;Fem.PotentialF];
   Fem.Log.Kin = [Fem.Kinetic;Fem.Kinetic];
   
   if norm(Fem.Utmp) == 0
         Fem.Log.Node{1} = Fem.Node0;
   else, Fem.Log.Node{1} = Fem.Node;
   end
   
   Fem.Log.Node{2} = Fem.Node;
   Fem.Log.U(1,:)  = Fem.Utmp.';
   Fem.Log.U(2,:)  = Fem.Utmp.';
   
   Fem.Log.Rotation{1} = R(:);
   Fem.Log.Rotation{2} = R(:);
   Fem.Log.Stretch{1}  = Q(:);
   Fem.Log.Stretch{2}  = Q(:);
   Fem.Log.Stress{1}   = Fem.VonMisesNodal;
   Fem.Log.Stress{2}   = Fem.VonMisesNodal;
   
   id = Fem.FindEdges('Hole');
   
   if ~isempty(id) && Fem.Dim == 2
       Vol = zeros(2,numel(id));
       for ii = 1:numel(id)
           Nds0 = Fem.Node0(id{ii},:);
           Nds = Fem.Node(id{ii},:);
           Vol(1,ii) = polyarea(Nds0(:,1),Nds0(:,2));
           Vol(2,ii) = polyarea(Nds(:,1),Nds(:,2));
       end
       Fem.Log.Volume = Vol;
   end
   
   if ~isempty(Fem.Flow)
       Fem.Log.z = [Fem.z0(:).';Fem.Log.AUX.z(:).'];
   end

   if ~isempty(Fem.Output)
       idNodes = Fem.Output(:,1);
       initVec = zeros(numel(idNodes),1);
       Fem.Log.Nx  = [initVec,Fem.Node(idNodes,1)];
       Fem.Log.Ny  = [initVec,Fem.Node(idNodes,2)];
       Fem.Log.Ux  = [initVec,ux(idNodes)];
       Fem.Log.Uy  = [initVec,uy(idNodes)];
       Fem.Log.Un  = [initVec,un(idNodes)];
       Fem.Log.Fx  = [initVec,fx(idNodes)];
       Fem.Log.Fy  = [initVec,fy(idNodes)];
       Fem.Log.Fn  = [initVec,fn(idNodes)];
       Fem.Log.Svm = [initVec,Fem.VonMisesNodal(idNodes)];
       Fem.Log.Sxx = [initVec,Fem.sxxNodal(idNodes)];
       Fem.Log.Syy = [initVec,Fem.syyNodal(idNodes)];
       Fem.Log.Sxy = [initVec,Fem.sxyNodal(idNodes)];
       Fem.Log.Exx = [initVec,Fem.exxNodal(idNodes)];
       Fem.Log.Eyy = [initVec,Fem.eyyNodal(idNodes)];
       Fem.Log.Exy = [initVec,Fem.exyNodal(idNodes)];
       Fem.Log.Exy = [initVec,Fem.exyNodal(idNodes)];
       Fem.Log.Rot = [{R(idNodes)}, {R(idNodes)}];
       Fem.Log.Str = [{Q(idNodes)}, {Q(idNodes)}];
   end
   
else
    
   Fem.Log.t   = vappend(Fem.Log.t,Fem.Time,1);
   Fem.Log.Psi = vappend(Fem.Log.Psi,Fem.Potential,1);
   Fem.Log.Vg  = vappend(Fem.Log.Vg ,Fem.PotentialG,1);
   Fem.Log.Vf  = vappend(Fem.Log.Vf ,Fem.PotentialF,1);
   Fem.Log.Kin = vappend(Fem.Log.Kin,Fem.Kinetic,1);
   Fem.Log.U   = vappend(Fem.Log.U,Fem.Utmp.',1);
   
   Fem.Log.Node{end + 1}     = Fem.Node;
   Fem.Log.Rotation{end + 1} = R(:);
   Fem.Log.Stretch{end + 1}  = Q(:);
   Fem.Log.Stress{end + 1}   = Fem.VonMisesNodal;
   
   
   id = Fem.FindEdges('Hole');
   if ~isempty(id) && Fem.Dim == 2
       Vol = zeros(1,numel(id));
       for ii = 1:numel(id)
           Nds = Fem.Node(id{ii},:);
           Vol(1,ii) = polyarea(Nds(:,1),Nds(:,2));
       end
       Fem.Log.Volume = vappend(Fem.Log.Volume,Vol,1);
   end
   
   if ~isempty(Fem.Flow)
       Fem.Log.z = vappend(Fem.Log.z,Fem.Log.AUX.z(:).',1);
   end

   if ~isempty(Fem.Output)
       idNodes = Fem.Output(:,1);
       Fem.Log.Nx   = vappend(Fem.Log.Nx,Fem.Node(idNodes,1),2);
       Fem.Log.Ny   = vappend(Fem.Log.Ny,Fem.Node(idNodes,2),2);
       Fem.Log.Ux   = vappend(Fem.Log.Ux,ux(idNodes),2);
       Fem.Log.Uy   = vappend(Fem.Log.Uy,uy(idNodes),2);
       Fem.Log.Un   = vappend(Fem.Log.Un,un(idNodes),2);
       Fem.Log.Fx   = vappend(Fem.Log.Fx,fx(idNodes),2);
       Fem.Log.Fy   = vappend(Fem.Log.Fy,fy(idNodes),2);
       Fem.Log.Fn   = vappend(Fem.Log.Fn,fn(idNodes),2);
       Fem.Log.Svm  = vappend(Fem.Log.Svm,Fem.VonMisesNodal(idNodes),2);
       Fem.Log.Sxx  = vappend(Fem.Log.Sxx,Fem.sxxNodal(idNodes),2);
       Fem.Log.Syy  = vappend(Fem.Log.Syy,Fem.syyNodal(idNodes),2);
       Fem.Log.Sxy  = vappend(Fem.Log.Sxy,Fem.sxyNodal(idNodes),2);
       Fem.Log.Exx  = vappend(Fem.Log.Exx,Fem.exxNodal(idNodes),2);
       Fem.Log.Eyy  = vappend(Fem.Log.Eyy,Fem.eyyNodal(idNodes),2);
       Fem.Log.Exy  = vappend(Fem.Log.Exy,Fem.exyNodal(idNodes),2);
       Fem.Log.Rot = [Fem.Log.Rot, {R(idNodes)}];
       Fem.Log.Str = [Fem.Log.Str, {Q(idNodes)}];
   end
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
%----------------------------------------------------------- mesh smoothing
function [Node,Node0] = UpdateNode(Fem,U)
    Node0 = Fem.get('Node0'); 
    Ntmp = Node0;

    u = FieldOperator(Fem,U);

    Ntmp(:,1) = Node0(:,1) + u(:,1);
    Ntmp(:,2) = Node0(:,2) + u(:,2);
    
    if Fem.Dim == 3
        Ntmp(:,3) = Node0(:,3) + u(:,3); 
    end

    Node = Ntmp;

    function u = FieldOperator(Fem,U)
        u  = zeros(Fem.NNode,Fem.Dim);
        
        for node = 1:Fem.NNode
            u(node,1) = U(Fem.Dim*node + (1-Fem.Dim),1);
            u(node,2) = U(Fem.Dim*node + (2-Fem.Dim),1);
            
            if Fem.Dim == 3
                u(node,3) = U(Fem.Dim*node,1);
            end
        end
    end
end
%----------------------------------------------------------- mesh smoothing
function Fem = UpdateSE3(Fem,U)
    
list = 1:Fem.NElem;
F2V  = Fem.Mesh.FaceToNode;

% for eacth element 
for ii = 1:Fem.NNode
    Rtmp = 0;
    Qtmp = 0;
    ElemFromNode = list(F2V(:,ii) > 0);
    
    for jj = ElemFromNode
        Qtmp = Qtmp + Fem.ElemStr{jj};
        Rtmp = Rtmp + Fem.ElemRot{jj};
    end
    
    [Ur,~,Vr] = svd(Rtmp);
    R{ii} = (Ur*Vr.');
    Q{ii} = Qtmp/numel(ElemFromNode);
end

Fem.Node     = UpdateNode(Fem,U);
Fem.Rotation = R;
Fem.Stretch  = R;

end
%---------------------------------------------- compute Hessian approximate
function z = UpdateAuxiliaryFlow(Fem)
z = Fem.Log.AUX.z;

% first EL-diff eval
f1 = Fem.Flow(Fem); 

Fem.Log.AUX.z = z + (2/3)*Fem.TimeStep*f1;
f2 = Fem.Flow(Fem);

% update integrands
z = z + 0.25*Fem.TimeStep*(f1 + 3*f2);
end

%%/////////////////////////////////////////////////// TOPOLOGY OPTIMIZATION
%------------------------------------------------------ objective function
function [f,dfdE,dfdV] = ObjectiveFunction(Fem)
    
u(:,1) = Fem.Utmp;
fDof = GetFreeDofs(Fem);
E = MaterialField(Fem);

if strcmp(Fem.OptimizationProblem,'Compliance') &&  ~Fem.Nonlinear
    f = (Fem.fInput).'*u;
    temp = cumsum(-u(Fem.i).*Fem.k.*u(Fem.j));
    temp = temp(cumsum(Fem.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliance') &&  Fem.Nonlinear
    K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
    f = (Fem.fExternal+Fem.fInput).'*u;
    lam = 0*u(:,1);
    lam(fDof) = K(fDof,fDof)\Fem.fInput(fDof);
    temp = cumsum(-u(Fem.i).*Fem.k.*lam(Fem.j));
    temp = temp(cumsum(Fem.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliant') && ~Fem.Nonlinear
    E = MaterialField(Fem);
    f = -Fem.OutputVector.'*u(:,1);
    K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
    lam = 0*u(:,1);
    lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
    temp = cumsum(u(Fem.i,1).*Fem.k.*lam(Fem.j,1));
    temp = temp(cumsum(Fem.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliant') && Fem.Nonlinear
    E = MaterialField(Fem);
    f = -Fem.OutputVector.'*u(:,1);
    K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
    %K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
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

if isempty(Fem.fnorm)
    Fem.fnorm = abs(norm(dfdz)); 
end

xval   = Fem.Density;
xmin   = zeros(N,1);
xmax   = ones(N,1);
f0val  = f;%(f/norm(f));
df0dx  = dfdz;%(dfdz/norm(dfdz));
df0dx2 = 0*df0dx;
fval   = g;
dfdx   = dgdz;
dfdx2  = dgdz*0;

A0 = 1; A = 0; C = 10000*ones(M,1); D = 0;

if iter == 1, Fem.low = xmin; Fem.upp = xmax; end
if iter > 1, Fem.xold1 = Fem.xold1; else, Fem.xold1 = 0; end
if iter > 2, Fem.xold2 = Fem.xold2; else, Fem.xold2 = 0; end

[xmma,~,~,~,~,~,~,~,~,Fem.low,Fem.upp] = mmasub(M,N,iter,xval,xmin,...
    xmax,Fem.xold1,Fem.xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,Fem.low,...
    Fem.upp,A0,A,C,D);

%if strcmp(Fem.OptimizationProblem,'Compliance')
%dx = clamp(xmma - xval,-Fem.ChangeMax*15,Fem.ChangeMax*15);
%else
dx = clamp(xmma - xval,-Fem.ChangeMax*15,Fem.ChangeMax*15);    
%end
zNew = xval + dx;

if iter >= 1, Fem.xold1 = xval; end
if iter >= 2, Fem.xold2 = Fem.xold1; end

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

if ~isempty(d)
    P = sparse(d(:,1),d(:,2),1-d(:,3)/R);
    P = spdiags(1./sum(P,2),0,Fem.NElem,Fem.NElem)*P;
else
    P = 1;
end
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
B5  = B(5); B6 = B(6); 
B56 = lerp(B5,B6,0.5);
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

sp1 = repmat(' ',[1,5-length(num2str(Fem.IterationMMA))]);
sp2 = repmat(' ',[1,5-length(num2str(Fem.Increment))]);
sp3 = repmat(' ',[1,6-length(num2str(Fem.Iteration))]);

fprintf([' %i',sp1,'| %i',sp2,'| %i',sp3,'| %1.3e | %0.3f | %1.3e | %0.3f |\n'],...
  Fem.IterationMMA,Fem.Increment,Fem.Iteration,abs(Fem.Objective(end)),...
  Fem.Constraint(end)+1,norm(Fem.Residual(FreeDofs)),norm(Fem.Change));

else
if Fem.Iteration == 1 && Fem.Increment == 1
fprintf(' Inc | Iter  | Residual  | Max. Svm  | Time | Beta  | dt     |\n');
fprintf('--------------------------------------------------------------\n');
end

sp1 = repmat(' ',[1,4-length(num2str(Fem.Increment))]);
sp2 = repmat(' ',[1,6-length(num2str(Fem.Iteration))]);

fprintf([' %i',sp1,'| %i',sp2,'| %1.3e | %1.3e | %1.2f | %1.3f | %1.3f  |\n'],...
    Fem.Increment,Fem.Iteration,norm(Fem.Residual(FreeDofs)),...
    max(Fem.s(:,1)),Fem.Time,Fem.LoadingFactor,Fem.TimeStep);   

end 

end
%------------------------------------------------------ optimization solver
function showInformation(Fem,Request)
if nargin < 2, Request = ''; end
    
if ~Fem.InformationBoolean
NodeNum = Fem.NElem;
NDofs   = numel(GetFreeDofs(Fem));
MaxNVer = max(cellfun(@numel,Fem.Element));      
MinNVer = min(cellfun(@numel,Fem.Element));      

switch(Request)
case('NonlinearFem')
fprintf('--------------------------------------------------------------\n');  
if NDofs>1e3, fprintf('* Free DOFs: %1.2fk \n',NDofs/1e3);
else, fprintf('* Free DOFs: %i \n',NDofs); end
if Fem.NElem>1e3, fprintf('* Elements: %1.2fk \n',NodeNum/1e3);
else, fprintf('* Element = %i \n',NodeNum); end
fprintf('* Element degree = P%1.0f-P%1.0f \n',MinNVer,MaxNVer);
fprintf('* Max iteration = %i \n', Fem.MaxIteration);
if Fem.Nonlinear, fprintf('* Nonlinear geometric = true\n');
else, fprintf('* Nonlinear geometric = false \n'); end
fprintf('* Solver time horizon = %i \n', Fem.TimeEnd);
fprintf('* Solver time step    = %i1.1e \n', Fem.TimeStep);
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
fprintf('\tRho  = %1.1e kg/mm2 \n', Fem.Material.Rho);
fprintf('\tZeta = %1.1e (-) \n', Fem.Material.Zeta);
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
            background(metropolis);
            if ~isempty(Fem.MovieAxis), axis(Fem.MovieAxis); end
            if ~isempty(Fem.MovieCAxis), caxis(Fem.MovieCAxis); end
            drawnow;
            gif(char(filename),'frame',gcf,'nodither');
            for ii = 1:10, gif; end
            
        otherwise
            background(metropolis);
            if ~isempty(Fem.MovieAxis), axis(Fem.MovieAxis); end
            if ~isempty(Fem.MovieCAxis), caxis(Fem.MovieCAxis); end
            drawnow;
            gif;
            if abs(Fem.Time-Fem.TimeEnd) <= 1e-3
               for ii = 1:14, gif; end
            end
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
%--------------------------------------------------- Kelvin-voight notation
function Sv = VoightNotation(S)
Sv = [S(1,1); S(2,2); S(3,3); S(1,2); S(2,3); S(1,3)]; 
end
%----------------------------------------------- compute von-Mises stresses
function [Svm, Svmm] = VonMises(S11,S22,S33,S12,S23,S13)
s11 = S11; s22 = S22; s33 = S33; s12 = S12; s23 = S23; s13 = S13;
Svm = sqrt(0.5*((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2 ...
    + 6*(s12.^2 + s23.^2 + s13.^2)));
Svmm = mean(Svm); 
end
%--------------------------------------------------- update centers of MESH
function Center = UpdateCenter(Fem)
    Center =  Fem.get('Center0'); 

    for el = 1:Fem.NElem
        E = Fem.Element{el};
        Center(el,:) = mean(Fem.Node(E,:),1);
    end

end
