classdef Fem < handle

    properties (Access = public)
        Mesh;
        Material;
        Density;
        Node;
        Element;
        NNode;
        NElem;
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
        MaxIterationMMA = 60;
        SolverStart = false;
        SolverStartMMA = false;
        SolverPlot = true;
        Assemble = false;
        Convergence = false; 
        Nonlinear = true;
        AssembledSystem = false;
        PrescribedDisplacement = false;
        VolumetricPressure = false;
        PressureLoad = 0;
        Type = 'PlaneStrain'
        LineStyle = '-';
        CMap = turbo;
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
        ChangeMax = 0.01;
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
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    obj = ConvertMeshToFem(obj,Mesh);
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
colormap(Fem.CMap);

switch(Request)
    case('0'),  Z = zeros(Fem.NNode,1);
    case('SDF'), Z = Fem.Mesh.SDF(Fem.Mesh.Node); Z = Z(:,end);
    case('Svm'), Z = Fem.VonMisesNodal;
    case('Sxx'), Z = Fem.sxxNodal;
    case('Syy'), Z = Fem.syyNodal;
    case('Sxy'), Z = Fem.sxyNodal;
    case('Exx'), Z = Fem.exxNodal;
    case('Eyy'), Z = Fem.eyyNodal;
    case('Exy'), Z = Fem.exyNodal;
    case('Fx'), Z = Fem.fxNodal;
    case('Fy'), Z = Fem.fyNodal;
    case('Fin'), [~,~,Z] = DisplacementField(Fem,Fem.fInternal);
    case('Fex'), [~,~,Z] = DisplacementField(Fem,Fem.fExternal);
    case('Un'), [~,~,Z] = DisplacementField(Fem,Fem.Utmp);
    case('Ux'), [Z,~,~] = DisplacementField(Fem,Fem.Utmp);
    case('Uy'), [~,Z,~] = DisplacementField(Fem,Fem.Utmp);
    case('E'), [~,~,Z] = MaterialField(Fem); S = 'flat'; 
    V = Fem.Node0; colormap(noir(-1)); background('w');
    case('E+'), [~,~,Z] = MaterialField(Fem); S = 'flat'; 
    colormap(noir(-1)); background('w');
    otherwise; flag = 1; Z = 0;
end

if length(Z) == Fem.NNode, Z = Smoothing(Fem,Z,1); end

if flag == 0
cla; axis equal;     
axis off; hold on; h{3} = [];
if length(Z) ~= Fem.NNode, T=Fem.Mesh.get('NodeToFace'); Z=T*Z; end

h{1} = patch('Faces',Fem.Mesh.get('Boundary'),'Vertices',Fem.Node0,...
    'LineStyle','-','Linewidth',1,'EdgeColor','k');

h{2} = patch('Faces',Fem.Mesh.get('ElemMat'),'Vertices',V,...
    'FaceVertexCData',Z,'Facecolor',S,'LineStyle',Fem.LineStyle,...
    'Linewidth',1.0,'FaceAlpha',1.0,'EdgeColor','k');

h{3} = patch('Faces',Fem.Mesh.get('Boundary'),'Vertices',V,...
    'LineStyle','-','Linewidth',2,'EdgeColor','k');

if Fem.VolumetricPressure
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density); 
    if strcmp(Request,'E+'), Pc = computeCentroid(Fem);
    else, Pc = Fem.Mesh.get('Center'); 
    end
    %h{4} = plot(Pc(id,1),Pc(id,2),'o','Color','k');
    h{5} = plot(Pc(id,1),Pc(id,2),'.','Color','k');
end
end

if flag == 1
switch(Request)
case('Residual'), semilogy(Fem.SolverResidual,'linewidth',1);
case('Stress'), semilogy(Fem.SolverVonMises,'linewidth',1);
case('Displace'), semilogy(Fem.SolverDisplace,'linewidth',1);
case('Objective'), plot(Fem.Objective,'linewidth',1);
case('Constraint'), plot(Fem.Constraint,'linewidth',1);
otherwise, flag = 2;
end
end

if flag == 2
    switch(Request)
    case('BC')
    if ~isempty(Fem.Support), Supp = Fem.Support; else, Supp = []; end
    if ~isempty(Fem.Load), Forc = Fem.Load; else, Forc = []; end
    if ~isempty(Fem.Spring), Sprg = Fem.Spring(:,1); else, Sprg = []; end
    if ~isempty(Fem.Output), Out = Fem.Output; else, Out = []; end
    if ~isempty(Fem.Pressure), Pres = Fem.Pressure(:,1); else, Pres = []; end
    
    hold on;
    for ii = 1:size(Supp,1)
        id = Supp(ii,1);
        if (Supp(ii,2) && Supp(ii,3)), symbol = '\Delta';
        else, symbol = '\Delta'; end
        plot(V(id,1),V(id,2),'marker','d','markersize',5,'Color',col(1));
    end
    
    for ii = 1:size(Forc,1)
        id = Forc(ii,1);
        if ispos(Forc(ii,2)), symbol = '\bf \rightarrow'; end
        if isneg(Forc(ii,2)), symbol = '\bf \leftarrow'; end
        if ispos(Forc(ii,3)), symbol = '\bf \uparrow'; end
        if isneg(Forc(ii,3)), symbol = '\bf \downarrow'; end
        vecmag = 0.15*sum(abs(axis))/4;
    end
    
    plotvector(V(Forc(:,1),:),Forc(:,2:3)*vecmag);
    
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
    if nargin < 3, ISOVALUE = 0.1; else, ISOVALUE = varargin{2};
    end
    showISO(Fem,ISOVALUE,0.5);
    end    
end

if Fem.Movie && strcmp(Request,'ISO')
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
   end
    
end

%-------------------------------------------------------------------- solve
function Fem = solve(Fem)
    
Fem.Log = [];    
Fem.TimeDelta = Fem.TimeEnd;
Fem.Increment = 1;
Fem.Iteration = 1;
Fem.Time = 0;
Fem.TimeDelta = 0;
Fem.EndIncrement = false;
Fem.Convergence = true;
Fem.TimeStep = Fem.TimeStep0 + 1e-6;
Fem.Density = clamp(Fem.Density,Fem.Ersatz,1);

Fem.Utmp = zeros(2*Fem.NNode,1);
Fem.Node = Fem.Node0;

showInformation(Fem,'NonlinearFem');

if Fem.SolverPlot || ~Fem.SolverStartMMA
    figure(101); Fem.show('Un');
end

while true 
    
    [Fem,Terminate] = SolverTimer(Fem);
    if Terminate, break; end
        
    flag = 0;
    Singular = false;
    Fem.Divergence = 0;
        
    % load increment
    Fem.LoadingFactor = sigmoid(Fem.Time/Fem.TimeEnd,0);
    
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
        Fem.SolverResidual = append(Fem.SolverResidual,...
            max(abs(Fem.Residual(FreeDofs))));
        
        Fem.SolverVonMises = append(Fem.SolverVonMises,...
            max(Fem.s(:,1)));
        
        Fem.SolverDisplace = append(Fem.SolverDisplace,...
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
            %[Pc,~,~,~] = ComputeCentroid(Fem);
            %Fem.Center = Pc;
%            Fem.Normal = Nv;
%            Fem.Edge = Ev;
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
           Fem.Log = cell(7,1);
           Fem.Log{1} = ux(idNodes); 
           Fem.Log{2} = uy(idNodes);
           Fem.Log{3} = un(idNodes);
           Fem.Log{4} = fx(idNodes); 
           Fem.Log{5} = fy(idNodes);
           Fem.Log{6} = fn(idNodes);
           Fem.Log{7} = Fem.VonMisesNodal(idNodes);
           Fem.Log{8} = Fem.sxxNodal(idNodes);
       else
           Fem.Log{1} = append(Fem.Log{1},ux(idNodes));
           Fem.Log{2} = append(Fem.Log{2},uy(idNodes));
           Fem.Log{3} = append(Fem.Log{3},un(idNodes));
           Fem.Log{4} = append(Fem.Log{4},fx(idNodes));
           Fem.Log{5} = append(Fem.Log{5},fy(idNodes));
           Fem.Log{6} = append(Fem.Log{6},fn(idNodes));
           Fem.Log{7} = append(Fem.Log{7},Fem.VonMisesNodal(idNodes));
           Fem.Log{8} = append(Fem.Log{8},Fem.sxxNodal(idNodes));
       end
               
    end
    
    if Fem.SolverPlot || ~Fem.SolverStartMMA
        figure(101); Fem.show('Svm'); drawnow;
    end
    
    if ~Fem.Nonlinear, break; end
end

end

%------------------------------------------------------ optimization solver
function Fem = optimize(Fem)
 
showInformation(Fem,'TopologyOptimization');
    
Fem.SpatialFilter = GenerateRadialFilter(Fem,Fem.FilterRadius);
Fem.SolverPlot = false;
Fem.IterationMMA = 0;
Fem.SolverStartMMA = true;
flag = true;
Visual = 'E+';

Fem.show(Visual);
drawnow;

while flag

    % update increment
    Fem.IterationMMA = Fem.IterationMMA + 1;
    
    % update material field
    [~,dEdy,~,dVdy] = MaterialField(Fem);
     
    % set load factor
    Fem.OptFactor = Fem.IterationMMA/Fem.MaxIterationMMA;
   
    % solve nonlinear finite elements
    Fem = Fem.solve();
    
    % objective and constraints
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
    Fem.Change = clamp(ZNew - Fem.Density,-Fem.ChangeMax,Fem.ChangeMax);
    Fem.Density =  Fem.Density + Fem.Change;
    
    % evaluate fitness
    Fem.Objective = append(Fem.Objective,f);       
    Fem.Constraint = append(Fem.Constraint,g);
    
    % check convergence
    [flag,Fem] = CheckConvergenceOpt(Fem);
    
    % draw
    Fem.show(Visual);
    drawnow;
end

Fem.SolverStartMMA = false;

end

%----------------------------------------------------- reconstruct topology
function Fem = former(Fem, Thickness)
Res = 200;    
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
SDF0 = GaussianFilter(SDF0,10);
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

%----------------------------------------------------- reconstruct topology
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
end
end

cla;
Uxx = scaleX*Rpt(1)*[min(min(min(X))) max(max(max(X)))];
Uyy = scaleY*Rpt(2)*[min(min(min(Y))) max(max(max(Y)))];
I = GaussianFilter((SDF >= varargin{1})*255,3);
[XX,YY] = meshgrid(linspace(Uxx(1),Uxx(2),size(I,1)),...
                   linspace(Uyy(1),Uyy(2),size(I,2)));

im = I(:); 
Ipos = reshape(im>0,size(I));
Ineg = reshape(im<=0,size(I));

[Dpos,~] = bwdist(Ipos,'euclidean');
[Dneg,~] = bwdist(Ineg,'euclidean');
Dist = rescale(Dpos) - rescale(Dneg);

%contourf(XX,YY,-Dist,[-1e-6,1e-6]);

image((max(Uxx)/max(Uxx))*rescale(Uxx),(max(Uyy)/max(Uxx))*rescale(Uyy),I);

axis equal; axis off; 
colormap(noir(-1)); 
caxis([0 1]);
background('w');
end

%----------------------------------------------------- reconstruct topology
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

%----------------------------------------------------------- toplogy design
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

%--------------------------------------------------------------- find nodes
function ElementList = FindElements(Fem,varargin)
    ElementList = FindNode(Fem.Mesh.Center,varargin{1:end});
end

%---------------------------------------------------------- add constraints
function Fem = AddConstraint(Fem,varargin)
    
for ii = 1:3:length(varargin)
  if size(varargin{ii+2},2) == 2
      if strcmp(varargin{ii},'PressureCell')
         Fem.VolumetricPressure = true; 
         Fem.(varargin{ii}) = [varargin{ii+1},repmat(varargin{ii+2},...
          [length(varargin{ii+1}),1])];
      elseif strcmp(varargin{ii},'Contact')
          Fem.(varargin{ii}) = {varargin{ii+1},varargin{ii+2}};
      else 
          BC = [varargin{ii+1},repmat(varargin{ii+2},...
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
function Fem = ConvertMeshToFem(Fem,Mesh)
Fem.Mesh     = Mesh;
Fem.Node     = Mesh.get('Node');
Fem.Node0    = Mesh.get('Node');
Fem.NNode    = Mesh.get('NNode');
Fem.Element  = Mesh.get('Element');
Fem.NElem    = Mesh.get('NElem');
Fem.Center   = Mesh.get('Center');
Fem.BdBox    = Mesh.get('BdBox');
Fem.Density  = ones(Fem.NElem,1);
Fem.Residual = zeros(2*Fem.NNode,1);
Fem.Utmp     = zeros(2*Fem.NNode,1);
Fem.TimeStep0 = Fem.TimeStep;

Fem.SpatialFilter = GenerateRadialFilter(Fem,Fem.FilterRadius);
[~,~,Fem.Normal,Fem.Edge] = computeCentroid(Fem);

Fem.ElemNDof = 2*cellfun(@length,Fem.Element);

if (~Fem.AssembledSystem && ~Fem.Nonlinear) || Fem.Nonlinear
Fem.i = zeros(sum(Fem.ElemNDof.^2),1);
Fem.j  = Fem.i; 
Fem.e  = Fem.i; 
Fem.fi = Fem.i; 
Fem.k  = Fem.i; 
Fem.m  = Fem.i; 
Fem.c  = Fem.i; 
Fem.t  = Fem.i; 
Fem.fb = Fem.i; 
Fem.ft = Fem.i;
Fem.s  = zeros(Fem.NNode,6); 
Fem.s  = zeros(Fem.NNode,3);
Fem.l  = zeros(Fem.NNode,1); 
Fem.v  = Fem.l;
end

end

%------------------------------------ assemble global finite-element system 
function Fem = AssembleGlobalSystem(Fem,ForceBuild)
    
if nargin < 2, ForceBuild = false; end

if (Fem.Iteration == 1 && Fem.Increment == 1)
    tmp = struct; tmp.Element = Fem.Element;
    tmp = TabulateShapeFunctions(tmp);
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
    eDof = reshape([2*Fem.Element{el}-1;
                    2*Fem.Element{el}],NDof,1);
    
    [Fe,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te] = ...
    Locals(Fem,Fem.Element{el},eDof,dV,V(el));
    
    ind1 = index+1:index+NDof^2;
    ind2 = index+1:index+NDof;
    ind3 = subindex+1:subindex+NDof/2;

    I = repmat(eDof ,1,NDof); J = I';
    Fem.e(ind1) = el;
    Fem.i(ind1) = I(:);
    Fem.j(ind1) = J(:);
    Fem.m(ind1) = Me(:);
    Fem.c(ind1) = Ce(:);
    Fem.k(ind1) = Ke(:);
    Fem.t(ind1) = Kte(:);
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
    
    index = index + NDof^2;
    subindex = subindex + NDof/2;
end
end

Fem.AssembledSystem = true;

end

%------------------------ assemble prescribed boundary forces/displacements
function Fem = AssembleBoundaryConditions(Fem)
    
[E,~,~,~] = MaterialField(Fem);

NSpring = size(Fem.Spring,1);
sp = sparse(2*Fem.NNode,2);
if NSpring ~=0
    sp(2*Fem.Spring(1:NSpring)-1) = Fem.Spring(1:NSpring,2);
    sp(2*Fem.Spring(1:NSpring)) = Fem.Spring(1:NSpring,3);
end
spMat = spdiags(sp(:),0,2*Fem.NNode,2*Fem.NNode);

NOutput = size(Fem.Output,1);
L = sparse(2*Fem.NNode,1);
if NOutput ~=0
    L(2*Fem.Output(1:NOutput,1)-1) = Fem.Output(1:NOutput,2);
    L(2*Fem.Output(1:NOutput,1)) = Fem.Output(1:NOutput,3);
end
    
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
Ktr = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
Fe = sparse(Fem.i,1,E(Fem.e).*Fem.fi);
F = sparse(2*Fem.NNode,1);

if Fem.Nonlinear, beta = Fem.OptFactor*Fem.LoadingFactor;
else, beta = 1;
end

if ~isempty(Fem.Load) && ~Fem.PrescribedDisplacement
NLoad = size(Fem.Load,1);
A = sign(Fem.Mesh.get('NodeToFace'));
Gamma = A(Fem.Load(1:NLoad,1),:);
W = sum(Gamma,2)/sum(sum(Gamma,2));
F(2*Fem.Load(1:NLoad,1)-1,1) = beta*W.*Fem.Load(1:NLoad,2);
F(2*Fem.Load(1:NLoad,1),1) = beta*W.*Fem.Load(1:NLoad,3);
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
    F(2*I(1:size(Y(I,:),1),1),1) = Uy;
    
    pDof = zeros(2*length(I),1);
    fDof = zeros(2*length(I),1);
    for kk = 1:length(I)
        pDof(2*kk-1) = 2*I(kk)-1; pDof(2*kk) = 2*I(kk);
        fDof(2*kk-1) = Ux(kk); fDof(2*kk) = Uy(kk);
    end
    
    Fem.PrescribedDisplacement = true;
end

if Fem.PrescribedDisplacement
    
    if ~isempty(Fem.Load)
    if abs(Fem.Load(1,2))>0
        pDof = reshape(2*Fem.Load(:,1)-1,NLoad,1);
        F(pDof) = Fem.Load(1:NLoad,2); id = 1;
    else 
        pDof = reshape(2*Fem.Load(:,1),NLoad,1);
        F(pDof) = Fem.Load(1:NLoad,3); id = 2;
    end
    end
    
    EMod = Fem.Material.Emod();
   
    if Fem.Nonlinear
        KTtmp = full(Ktr);
        I = eye(2*Fem.NNode,2*Fem.NNode);
        KTtmp(pDof,:) = 0*I(pDof,:);
        KTtmp(pDof,pDof) = -EMod*eye(length(pDof));
        Ktr = KTtmp;

        if isempty(Fem.Contact), alpha = Fem.TimeStep*Fem.OptFactor;
        else, alpha = 1; end
        
        if Fem.Iteration == 1
        F(pDof) = Fe(pDof)-EMod*F(pDof)*alpha;
        else, F(pDof) = Fe(pDof); 
        end
    else
        Ktmp = full(K); I = eye(2*Fem.NNode,2*Fem.NNode);
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
            %ntmp = w(:)'*N{jj};
            ee = Fem.Edge{jj};
            etmp(kk) = max(ee(Elem{jj} == IdNode));
            kk = kk+1;
        end
        
        ntmp = 0.5*(ntmp/length(el));
        
        PressureForce(ii,:) = ntmp*Pressmag;
    end
    
    NPLoad = size(Fem.Pressure,1);
    F(2*Fem.Pressure(1:NPLoad,1)-1,1) = beta*PressureForce(1:NPLoad,1)/NPLoad;
    F(2*Fem.Pressure(1:NPLoad,1),1) = beta*PressureForce(1:NPLoad,2)/NPLoad;
end

if Fem.VolumetricPressure
    
    Pc = Fem.Mesh.get('Center');
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
    idnull = setdiff(1:Fem.NElem,id);
    %[~,~,V,~] = MaterialField(Fem);
    %[~,~,W,~] = MaterialField(Fem,true);
    W = ones(Fem.NElem,1);
    W(idnull) = 0;
    %W = Fem.SpatialFilter*W;
    %V = Fem.Mesh.get('NodeToFace')*V;
    %Ev = wthresh(kron(V,[1;1]),'h',-Fem.VoidTolerance);
    %Ev = wthresh(V,'h',Fem.VoidTolerance);
    %W(idnull) = 0;
    %W = Fem.SpatialFilter*W;
    %fe = beta*Ev.*sparse(Fem.i,1,W(Fem.e).*Fem.ft);
    %W = Fem.SpatialFilter*W;
    %E = wthresh(Fem.Mesh.get('NodeToFace').'*E,'s',Fem.VoidTolerance);
    %V = Fem.Mesh.get('NodeToFace')*V;
    %Ev = wthresh(kron(V,[1;1]),'h',Fem.VoidTolerance);
    
    %Ft = beta*sparse(Fem.i,1,W(Fem.e).*V(Fem.e).*Fem.ft);
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
            %[W,~,V,~] = MaterialField(Fem);
            %V = Fem.Mesh.get('NodeToFace')*V;
            %Ev = wthresh(kron(V,[1;1]),'h',Fem.VoidTolerance);
            W = ones(Fem.NElem,1);
            W(idnull) = 0;
            %W = Fem.SpatialFilter*W;
            %fe = beta*Ev.*sparse(Fem.i,1,W(Fem.e).*Fem.ft);
            fe = sparse(Fem.i,1,W(Fem.e).*Fem.ft);
            Fem.dFdE(:,ii) = ((fe - Ft)/(dz));
        end
        Fem.Density = z0;
    end
else
    Ft = sparse(Fem.NNode*2,1);
end

Fem.fInternal = Fe; 
Fem.fExternal = F - spMat*Fem.Utmp + Ft;
Fem.TangentStiffness = Ktr + spMat;
Fem.Stiffness = K + spMat;
Fem.Residual = Fem.fInternal - Fem.fExternal;
Fem.OutputVector = L;

end

%---------------------------------------------- get free degrees-of-freedom
function FreeDofs = GetFreeDofs(Fem)
NSupp = size(Fem.Support,1);
FixedDofs = [Fem.Support(1:NSupp,2).*(2*Fem.Support(1:NSupp,1)-1);
    Fem.Support(1:NSupp,3).*(2*Fem.Support(1:NSupp,1))];
FixedDofs = FixedDofs(FixedDofs>0);
AllDofs = 1:2*Fem.NNode;

FreeDofs = setdiff(AllDofs,FixedDofs);
end

%--------------------------------------------------- local element matrices
function [Fe,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te] = Locals(Fem,eNode,eDof,dV,Rb)
% get order
nn  = length(eNode);
    
% get gauss points and weights
W  = Fem.ShapeFnc{nn}.W;
Q  = Fem.ShapeFnc{nn}.Q;    
Fe  = zeros(2*nn,1);
Te  = zeros(2*nn,1);
Me  = zeros(2*nn,2*nn);
Ce  = zeros(2*nn,2*nn);
Ke  = zeros(2*nn,2*nn);
Kte = zeros(2*nn,2*nn);
SGP = zeros(length(W),6);
EGP = zeros(length(W),6);
Qe  = ones(nn,1);

Nshp = length(Fem.ShapeFnc{nn}.N(:,:,1));
NNe = zeros(length(W),Nshp);

Et = [dV;dV;0];

% quadrature loop
for q = 1:length(W)
    dNdxi = Fem.ShapeFnc{nn}.dNdxi(:,:,q);
    N = Fem.ShapeFnc{nn}.N(:,:,q);
    J0 = Fem.Node0(eNode,:).'*dNdxi;
    dNdx = dNdxi/J0;
    dJ = (det(J0));
       
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
    
    if strcmp(Fem.Type,'PlaneStress')
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

%Ke = 0.5*(Ke+Ke');
%Kte = 0.5*(Kte+Kte');
end

%----------------------------------------------- TRIANGULAR SHAPE FUNCTIONS
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

%----------------------------------------------- TRIANGULAR SHAPE FUNCTIONS
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
function F = DeformationGradient(~,U,dNdx)
nn = length(U)/2;
UU = zeros(nn,2);
UU(:,1) = U(1:2:2*nn);
UU(:,2) = U(2:2:2*nn);
F = (dNdx'*UU)';
F = [F(1,1)+1,F(1,2),0; F(2,1),F(2,2)+1,0;0,0,1];
end

%---------------------------------------------------------- polar decompose
function [R,S,V] = PolarDecomposition(~,F)
[M, N] = size(F);
if M ~= N
    error('Matrix must be square.');
end
C = F'*F;
[Q0, lambdasquare] = eig(C);
lambda = sqrt(diag((lambdasquare))); % extract the components
% Uinv is the inverse of U and is constructed with the help of Q0. Uinv is
% produced in the same base as F not in the base of its eigenvectors.
Uinv = repmat(1./lambda',size(F,1),1).*Q0*Q0';
% Using the definition, R, U and V can now be calculated
R = F*Uinv;
S = R'*F;
V = F*R';
end

%---------------------------------------------------------- strain operator
function [Bn,Bg,NN,tau] = NonlinearStrainOperator(~,N,dNdx,F)
nn = length(N);

NN = zeros(2,2*nn);
NN(1,1:2:2*nn) = N.';
NN(2,2:2:2*nn) = N.';
dNdxX = dNdx(:,1).';
dNdxY = dNdx(:,2).';

%if strcmp(Fem.Type,'PlaneStress') || strcmp(Fem.Type,'PlaneStrain')
    Bn = zeros(3,2*nn); 
    Bg = zeros(4,2*nn);
    
    Bn(1,1:2:2*nn) = dNdxX*F(1,1);
    Bn(1,2:2:2*nn) = dNdxX*F(2,1);
    Bn(2,1:2:2*nn) = dNdxY*F(1,2);
    Bn(2,2:2:2*nn) = dNdxY*F(2,2);
    Bn(3,1:2:2*nn) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(3,2:2:2*nn) = dNdxX*F(2,2) + dNdxY*F(2,1);
   
    Bg(1,1:2:2*nn) = dNdxX;
    Bg(2,1:2:2*nn) = dNdxY;
    Bg(3,2:2:2*nn) = dNdxX;
    Bg(4,2:2:2*nn) = dNdxY;
    
    tau = 1;
%end

end 

%------------------------------------------------------ isotropic reduction
function [S, D, G] = IsotropicReduction(~,D0,S0)

D = zeros(3,3); 
S = zeros(3,1);

% if Type == 1
% D(1:2,1:2) = D0(1:2,1:2); D(3,3) = D0(4,4);
% SIG = [S0(1),S0(4);
%        S0(4),S0(2)];
% G = kron(eye(2),SIG);
% S(1) = S0(1); 
% S(2) = S0(2); 
% S(3) = S0(4); 
% 
% elseif Type == 2
D(1:2,1:2) = D0(1:2,1:2); D(3,3) = D0(4,4);
SIG = [S0(1),S0(4);
       S0(4),S0(2)];
G = kron(eye(2),SIG);
S(1) = S0(1); 
S(2) = S0(2); 
S(3) = S0(4); 

% elseif strcmp(Fem.Type,'AxiSymmetric')
% 
% end
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
    if round((100*Fem.Iteration/Fem.MaxIteration)) > 75
        Fem.MaxIteration = round(Fem.MaxIteration)*1.25;
    end
end

if ~Fem.SolverStartMMA, ProcessMonitor(Fem); end

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
    Fem.MaxIteration = Fem.MaxIteration + 15;
    if ~Fem.SolverStartMMA
    fprintf('------------------------------------------------------------|\n');
    fprintf(' Bisection                                                  |\n');
    fprintf('------------------------------------------------------------|\n');
    end
end

if ~Fem.SolverStart, Fem.Increment = 1; Fem.SolverStart = true;
else, Fem.Increment = Fem.Increment + 1;
end

if (Fem.SolverStartMMA && Terminate)||(Fem.SolverStartMMA && ~Fem.Nonlinear)
    ProcessMonitor(Fem); 
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
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
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
%f0val = (f/Fem.fnorm);
%df0dx = (dfdz/Fem.fnorm);
%f0val = f;
%df0dx = dfdz;
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

if strcmp(Fem.OptimizationProblem,'Compliance'), zNew = xmma;
else
alpha = max(0.9-iter/Fem.MaxIterationMMA,0);
zNew = (1-alpha)*xmma + alpha*xval;
end

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

for el = 1:size(PS1,1)   
    dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    
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
if Fem.IterationMMA == 1 && Fem.Iteration == 1
fprintf(' Iter | Inc  | Residual  | Obj. fun. | Vol.  | Change | p    |\n');
fprintf('--------------------------------------------------------------\n');
end

fprintf(' %1.0f\t  | %1.0f\t | %1.3e | %1.3e | %0.3f | %0.3f  | %0.2f \n',...
    Fem.IterationMMA,Fem.Increment,norm(Fem.Residual(FreeDofs)),...
    (Fem.Objective(end)),Fem.Constraint(end)+1,norm(Fem.Change),Fem.Penal);

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
fprintf('==============================================================\n');    
if Fem.NElem>1e3, fprintf('Number of elem: %1.2fk \n',NodeNum/1e3);
else, fprintf('Number of elem: %i \n',NodeNum); end
fprintf('Element type: P%1.0f-P%1.0f \n',MinNVer,MaxNVer);
fprintf('Max. iterations: %i \n', Fem.MaxIteration);
if Fem.Nonlinear, fprintf('Nonlinear geom: true\n');
else, fprintf('Nonlinear geom: false \n'); end
showMaterialInfo(Fem)
fprintf('==============================================================\n');
    
case('TopologyOptimization')
fprintf('==============================================================\n');    
if Fem.NElem>1e3, fprintf('Number of elem: %1.2fk \n',NodeNum/1e3);
else, fprintf('Number of elem: %i \n',NodeNum); end
fprintf('Element type: P%1.0f-P%1.0f \n',MinNVer,MaxNVer);
fprintf('Filter radius: %1.2f \n', Fem.FilterRadius);
fprintf('Max. iterations: %i \n', Fem.MaxIterationMMA);
fprintf('Optimization: %s \n', Fem.OptimizationProblem);
fprintf('Material scheme: %s \n', Fem.MaterialInterpolation);
if Fem.Nonlinear, fprintf('Nonlinear geom: true\n');
else, fprintf('Nonlinear geom: false \n'); end
showMaterialInfo(Fem)
fprintf('==============================================================\n');
end

Fem.InformationBoolean = true;
end

end

%------------------------------------------------- show material properties
function showMaterialInfo(Fem)
if strcmp(Fem.Material.Type,'Mooney'), fprintf('Material: Mooney \n');
c1 = Fem.Material.C10; c2 = Fem.Material.C01; d = Fem.Material.K;
fprintf('C10 = %1.e, C01 = %1.3e, K = %1.3e \n', c1,c2,d);
elseif strcmp(Fem.Material.Type,'Yeoh'), fprintf('Material = Yeoh \n');
c1 = Fem.Material.C1*1e3; c2 = Fem.Material.C2*1e3; c3 = Fem.Material.C3*1e3;
fprintf('C1 = %1.1e kPa, C2 = %1.1e kPa, C3 = %1.1e kPa\n', c1,c2,c3);
elseif strcmp(Fem.Material.Type,'NeoHookean')
fprintf('Material = Neo-Hookean \n');
E = Fem.Material.E*1e3; Nu = Fem.Material.Nu; 
fprintf('E = %1.1e kPa, Nu = %1.2f [-] \n', E,Nu);
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
