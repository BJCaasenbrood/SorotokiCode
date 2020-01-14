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
    end
    
    properties (Access = private)
        Support;
        Load;
        Spring;
        Output;
        PressureCell;
        Contraction;
        FixedDensity;
        ElemNDof;
        Normal;
        VonMisesNodal;
        sxxNodal; 
        syyNodal;
        sxyNodal;
        fxNodal;
        fyNodal;
        SPxyNodal;
        ElemMat;
        Rho = 1e12;
        Zeta = 0.1;
        fInternal = rand(1)*1e-8;
        fExternal = rand(1)*1e-8;
        Stiffness;
        TangentStiffness;
        Residual = Inf;
        Time = 0;
        TimeEnd = 1.0;
        TimeStep = 0.1;
        TimeStep0 = 0.1;
        TimeDelta = 0;
        TimeStepMin = 1e-3;
        EndIncrement;
        LoadingFactor;
        ResidualNorm = 1e-4;
        IterationMMA = 1;
        Iteration = 1;
        Increment = 1;
        ShapeFnc;
        U = 0;
        Utmp = 0;
        Node0;
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
        LineStyle = 'none';
        I3 = eye(3); O3 = zeros(3);
        i; j; m; fi; t; e; c; s; v; l; k; fb; ed; fb0;
        VolumeInfill = 0.3;
        Ersatz = 1e-3; 
        Penal = 1;
        PenalMax = 5;
        PenalStep = 10;
        ChangeMax = 10;
        Beta = 1.5;
        FilterRadius = 0.5;
        SpatialFilter;
        OptFactor = 1;
        MaterialInterpolation = 'SIMP';
        OptimizationProblem = 'Compliance';
        xold1; xold2; upp; low;
        fnorm;
        zMin = 0;
        zMax = 1;
        OutputVector;
        Change;
        dFdE;
        Periodic;
               
        VoidTolerance = 0.05; 
        
        CutOff = 0.84;
        Threshold = .003;
        GaussianWeight = 10;
        GaussianRadius = .45;
        Reflect;
        NumGridX; NumGridY
        Grid; Grideval;
        Obj = rand(1)*1e-8;
        Con = rand(1)*1e-8;
        
        
        Colorbar = false;
        ColorbarLoc = 'southoutside';
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Fem(Mesh,varargin) 
    obj = ConvertMeshToFem(obj,Mesh);
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
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
    
if nargin<2, Request = -1; 
else, Request = varargin{1}; end

Shading = 'interp'; V = Fem.Node;

switch(Request)
    case('Svm'), Z = Fem.VonMisesNodal;
    case('Sxx'), Z = Fem.sxxNodal;
    case('Syy'), Z = Fem.syyNodal;
    case('Sxy'), Z = Fem.sxyNodal;
    case('Fx'), Z = Fem.fxNodal;
    case('Fy'), Z = Fem.fyNodal;
    case('Fi'), [~,~,Z] = DisplacementField(Fem,Fem.fInternal);
    case('Un'), [~,~,Z] = DisplacementField(Fem,Fem.Utmp);
    case('Ux'), [Z,~,~] = DisplacementField(Fem,Fem.Utmp);
    case('Uy'), [~,Z,~] = DisplacementField(Fem,Fem.Utmp);
    case('E'), [~,~,Z] = MaterialField(Fem); %Z = 1-Z.^0.75;
        Shading = 'flat'; V = Fem.Node0; colormap(bluesea(-1)); 
        figure(101); background('w');
    otherwise; Z = Fem.VonMisesNodal;
end

cla; axis equal; axis off; hold on; h{3} = [];

if length(Z) ~= Fem.NNode, T=Fem.Mesh.get('NodeToFace'); Z=T*Z; end

h{1} = patch('Faces',Fem.Mesh.get('Boundary'),'Vertices',Fem.Node0,...
    'LineStyle','-','Linewidth',1,'EdgeColor','k');

h{2} = patch('Faces',Fem.Mesh.get('ElemMat'),'Vertices',V,...
    'FaceVertexCData',Z,'Facecolor',Shading,'LineStyle',Fem.LineStyle,...
    'Linewidth',1.0,'FaceAlpha',1.0,'EdgeColor','k');

h{3} = patch('Faces',Fem.Mesh.get('Boundary'),'Vertices',V,...
    'LineStyle','-','Linewidth',2,'EdgeColor','k');

if Fem.VolumetricPressure
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density); 
    Pc = Fem.Mesh.get('Center');
    h{4} = plot(Pc(id,1),Pc(id,2),'o','Color',lightblue);
    h{5} = plot(Pc(id,1),Pc(id,2),'.','Color',lightblue);
end

end
%-------------------------------------------------------------------- solve
function showBC(Fem)
    
vert = Fem.Node;
Supp = Fem.Support(:,1);
if ~isempty(Fem.Load), Forc = Fem.Load(:,1); else, Forc = []; end
if ~isempty(Fem.Spring), Sprg = Fem.Spring(:,1); else, Sprg = []; end
if ~isempty(Fem.Output), Out = Fem.Output(:,1); else, Out = []; end

hold on;
for ii = 1:length(Supp)
    id = Supp(ii);
    plot(vert(id,1),vert(id,2),'+','markersize',10,'Color','r');
end

for ii = 1:length(Forc)
    id = Forc(ii);
    plot(vert(id,1),vert(id,2),'d','markersize',10,'Color','b');
end

for ii = 1:length(Out)
    id = Out(ii);
    plot(vert(id,1),vert(id,2),'d','markersize',10,'Color','b');
end

for ii = 1:length(Sprg)
    id = Sprg(ii);
    plot(vert(id,1),vert(id,2),'>','markersize',10,'Color','g');
end


end

%-------------------------------------------------------------------- solve
function Fem = solve(Fem)
    
Fem.TimeDelta = Fem.TimeEnd;
Fem.Increment = 1;
Fem.Iteration = 1;
Fem.Time = 0;
Fem.TimeDelta = 0;
Fem.EndIncrement = false;
Fem.Convergence = true;
Fem.Utmp = zeros(2*Fem.NNode,1);
Fem.TimeStep = Fem.TimeStep0 + 1e-6;
Fem.Node = Fem.Node0;
Fem.Utmp = zeros(2*Fem.NNode,1);

if ~Fem.SolverStartMMA, showInformation(Fem); end

while true
    
    [Fem,Terminate] = SolverTimer(Fem);
    if Terminate, break; end
        
    flag = 0;
    Singular = false;
    
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
        
        if (Fem.Iteration == 1 && Fem.Increment == 1) 
            Delta = zeros(2*Fem.NNode,1);
        else
            Delta = Fem.Utmp;
        end
           
        if rcond(full(A)) >= 1e-20, DeltaU = A\B;
        else, Singular = true; DeltaU = Fem.Utmp(FreeDofs)*0;
        end
            
        if Fem.Nonlinear, Delta(FreeDofs,1)=Delta(FreeDofs,1)-DeltaU(:,1);
        else, Delta(FreeDofs,1) = DeltaU(:,1); B = Fem.ResidualNorm; end
        
        % update node
        [Fem.Node,~] = UpdateNode(Fem,Delta);
        
        % check convergence
        [flag,Fem] = CheckConvergence(Fem,B,Singular); 
        
        Fem.Utmp = Delta;
        
        %Fem.Center = ComputeCentroid(Fem);
        %Fem.Normal = ComputeNormal(Fem);
        Fem.Iteration = Fem.Iteration + 1;

    end 
    
    if ~Fem.SolverStartMMA
        Fem.VonMisesNodal = full(sparse(Fem.l,1,Fem.s(:,1))./sparse(Fem.l,1,Fem.v));
        Fem.sxxNodal = full(sparse(Fem.l,1,Fem.s(:,2))./sparse(Fem.l,1,Fem.v));
        Fem.syyNodal = full(sparse(Fem.l,1,Fem.s(:,3))./sparse(Fem.l,1,Fem.v));
        Fem.sxyNodal = full(sparse(Fem.l,1,Fem.s(:,4))./sparse(Fem.l,1,Fem.v));
        force = full(sparse(Fem.i,1,Fem.fi));
        Fem.fxNodal = force(2*(1:Fem.NNode)-1);
        Fem.fyNodal = force(2*(1:Fem.NNode));
    end
    
    if Fem.SolverPlot || ~Fem.SolverStartMMA, figure(101); Fem.show('Svm'); end
    
    if ~Fem.Nonlinear, break; end
end

end

%------------------------------------------------------ optimization solver
function Fem = optimize(Fem)
 
showInformation(Fem);
    
Fem.SpatialFilter = GenerateRadialFilter(Fem,Fem.FilterRadius);
Fem.SolverPlot = false;
Fem.IterationMMA = 0;
Fem.SolverStartMMA = true;
flag = true;

fig = Fem.show('E');
filename = string(['topo_', char(datetime(now,'ConvertFrom','datenum')),'.gif']);
filename = erase(filename,[":"," "]);
gif(char(filename));

while flag

    % update increment
    Fem.IterationMMA = Fem.IterationMMA + 1;
    
    % update material field
    [~,dEdy,~,dVdy] = MaterialField(Fem);
     
    % solve nonlinear finite elements
    Fem.OptFactor = 1;
   
    % solve nonlinear finite elements
    Fem = Fem.solve();
    
    % compute cost functionals and analysis sensitivities
    [f,dfdE,dfdV] = ObjectiveFunction(Fem);
    [g,dgdE,dgdV] = ConstraintFunction(Fem);
    
    % compute design sensitivities
    dfdz = Fem.SpatialFilter'*(dEdy.*dfdE + dVdy.*dfdV);
    dgdz = Fem.SpatialFilter'*(dEdy.*dgdE + dVdy.*dgdV);
    
    % compute design variable
    
    %[Fem,ZNew] = UpdateSchemeMMA(Fem,f,dfdz,g,dgdz);
    [Fem,ZNew] = UpdateSchemeOC(Fem,dfdz,g,dgdz);
    
    % determine material change
    Fem.Change = clamp(ZNew - Fem.Density,-Fem.ChangeMax,Fem.ChangeMax);
    Fem.Density =  Fem.Density + Fem.Change;
    
    [flag,Fem] = CheckConvergenceOpt(Fem);
    
    Fem.Obj = f; Fem.Con = g;
    
    set(fig{2},'FaceVertexCData',Fem.SpatialFilter*Fem.Density); 
    
    if Fem.VolumetricPressure
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density);    
    Pc = Fem.Mesh.get('Center');   
    set(fig{4},'XData',Pc(id,1),'YData',Pc(id,2)); 
    set(fig{5},'XData',Pc(id,1),'YData',Pc(id,2)); 
    end
    
    drawnow;
    gif;
end

gif('close'); 
Fem.SolverStartMMA = false;

end

%----------------------------------------------------- reconstruct topology
function Fem = former(Fem)

Res = 200;    
    
verts = Fem.Node0; 
V = Fem.Mesh.get('NodeToFace')*Fem.SpatialFilter*Fem.Density;


x = verts(:,1); y = verts(:,2);
xq = linspace(Fem.BdBox(1),Fem.BdBox(2),Res); 
yq = linspace(Fem.BdBox(3),Fem.BdBox(4),Res);

[xq,yq] = meshgrid(xq,yq);

% [xq,yq] = ndgrid(xq,yq);
% F = griddedInterpolant(xq,yq,V);
% [a,b] = ndgrid(m,n);
% c = F(a,b);

P = griddata(x,y,V,xq,yq);
Dist = Fem.Mesh.SDF([xq(:),yq(:)]); Dist = Dist(:,end);
P = P(:); P(Dist > 1e-3) = 0; 
P = (reshape(P,Res,Res));

Fem.Topology = cat(3,xq,yq,P);

end

%----------------------------------------------------- reconstruct topology
function Fem = showTopo(Fem)
    
f = figure;

ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]); 

b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',0.5, 'min',Fem.Ersatz, 'max',1-Fem.Ersatz);
          
bgcolor = f.Color;

bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,30,100,20],...
                'String','Dilation factor','BackgroundColor',bgcolor,...
                'FontWeight','bold');
            
V = Fem.Topology;
            
b.Callback = @(es,ed) cb(es,ed,Fem);

    function cb(es,ed,Fem)
        vert = Fem.Node0;

%         contour(V(:,:,1),V(:,:,2),V(:,:,3),[es.Value es.Value],...
%              'linestyle','-','EdgeColor','k','Linewidth',2), 
        
        Z = GaussianFilter(V(:,:,3),5);

        contourf(V(:,:,1),V(:,:,2),Z,[es.Value es.Value],...
            'linestyle','-','EdgeColor','k','Linewidth',2), 
        
        patch('Faces',Fem.Mesh.get('Boundary'),'Vertices',vert,...
         'LineStyle','-','Linewidth',2,'EdgeColor','k');
     
        axis equal 
        axis off;
        colormap(inferno);
        caxis([0 1]);
    end

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
Fem.Density = Fem.SpatialFilter*(1.2*Fem.Density*Fem.VolumeInfill);
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
      Constraint = [varargin{ii+1},repmat(varargin{ii+2},...
          [length(varargin{ii+1}),1])];
      Fem.(varargin{ii}) = [Fem.(varargin{ii});Constraint];
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
Fem.Mesh = Mesh;
Fem.Node = Mesh.get('Node');
Fem.Node0 = Mesh.get('Node');
Fem.NNode = Mesh.get('NNode');
Fem.Element = Mesh.get('Element');
Fem.NElem = Mesh.get('NElem');
Fem.BdBox = Mesh.get('BdBox');
Fem.Density = ones(Fem.NElem,1);
Fem.SpatialFilter = GenerateRadialFilter(Fem,Fem.FilterRadius);
Fem.TimeStep0 = Fem.TimeStep;
Fem.Center = Fem.Mesh.get('Center');
Fem.Normal = ComputeNormal(Fem);
Fem.Residual = zeros(2*Fem.NNode,1);
end

%------------------------------------ assemble global finite-element system 
function Fem = AssembleGlobalSystem(Fem)

Fem.ElemNDof = 2*cellfun(@length,Fem.Element);
if (~Fem.AssembledSystem && ~Fem.Nonlinear) || Fem.Nonlinear
Fem.i = zeros(sum(Fem.ElemNDof.^2),1);
Fem.j = Fem.i; Fem.e = Fem.i; Fem.fi = Fem.i; Fem.k = Fem.i; 
Fem.m = Fem.i; Fem.c = Fem.i; Fem.t = Fem.i; Fem.fb = Fem.i; 
Fem.fb0 = Fem.i;
Fem.s = zeros(Fem.NNode,7);
Fem.l = zeros(Fem.NNode,1); Fem.v = Fem.l;
end

if (Fem.Iteration == 1 && Fem.Increment == 1)
tmp = struct; tmp.Element = Fem.Element;
tmp = TabulateShapeFunctions(tmp);  
Fem.ShapeFnc = tmp.ShapeFnc;
end

if Fem.VolumetricPressure
   [E,~,~,~] = MaterialField(Fem);
   E = Fem.Mesh.get('NodeToFace')*E;
   %dE = Fem.Mesh.get('NodeToFace')*dE;
   Ev = kron(E,[1;1]);
   %dEv = kron(dE,[1;1]);
end

index = 0; subindex = 0;

if (~Fem.AssembledSystem && ~Fem.Nonlinear) || Fem.Nonlinear
    
for el = 1:Fem.NElem
    
    NDof = Fem.ElemNDof(el);
    eDof = reshape([2*Fem.Element{el}-1;
        2*Fem.Element{el}],NDof,1);
    
    [Fe,Qe,~,~,Ke,Kte,Svme,SS] = ...
    Locals(Fem,Fem.Element{el},eDof);

    Ve = ElemCellVolumeForce(Fem,Fem.Element{el},el);   
    
    I = repmat(eDof ,1,NDof); J = I';
    Fem.e(index+1:index+NDof^2) = el;
    Fem.i(index+1:index+NDof^2) = I(:);
    Fem.j(index+1:index+NDof^2) = J(:);
    %Fem.m(index+1:index+NDof^2) = Me(:);
    %Fem.c(index+1:index+NDof^2) = Ce(:);
    Fem.k(index+1:index+NDof^2) = Ke(:);
    Fem.t(index+1:index+NDof^2) = Kte(:);
    Fem.fi(index+1:index+NDof) = Fe(:);
    Fem.fb0(index+1:index+NDof) = Kte*Ve(:);
    
    if Fem.VolumetricPressure
        Ed = diag(Ev(eDof));
        Fem.ed(index+1:index+NDof^2) = Ed(:);
        Fem.fb(index+1:index+NDof) = Ed*Ke*Ve(:);
    end
    
    Fem.s(subindex+1:subindex+NDof/2,1) = Svme(:);
    Fem.s(subindex+1:subindex+NDof/2,2) = SS(:,1);
    Fem.s(subindex+1:subindex+NDof/2,3) = SS(:,2);
    Fem.s(subindex+1:subindex+NDof/2,4) = SS(:,4);
    Fem.v(subindex+1:subindex+NDof/2) = Qe(:);
    Fem.l(subindex+1:subindex+NDof/2) = Fem.Element{el}(:);
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
F  = sparse(2*Fem.NNode,1);

if Fem.Nonlinear, beta = Fem.OptFactor*Fem.LoadingFactor;
else, beta = 1;
end

if ~isempty(Fem.Load)
NLoad = size(Fem.Load,1);
F(2*Fem.Load(1:NLoad,1)-1,1) = beta*Fem.Load(1:NLoad,2)/NLoad;
F(2*Fem.Load(1:NLoad,1),1) = beta*Fem.Load(1:NLoad,3)/NLoad;
end

if Fem.PrescribedDisplacement
    F = sparse(2*Fem.NNode,1);
    if abs(Fem.Load(1,2))>0
        pDof = reshape(2*Fem.Load(:,1)-1,NLoad,1);
        F(pDof) = Fem.Load(1:NLoad,2); id = 1;
    else
        pDof = reshape(2*Fem.Load(:,1),NLoad,1);
        F(pDof) = Fem.Load(1:NLoad,3); id = 2;
    end
    
    if strcmp(Fem.Material.Type,'Mooney'), EMod = Fem.Material.C10;
    elseif strcmp(Fem.Material.Type,'Yeoh'), EMod = 6*Fem.Material.C1;
    elseif strcmp(Fem.Material.Type,'Linear'), EMod = Fem.Material.E;
    end
   
    if Fem.Nonlinear
        KTtmp = full(Ktr);
        I = eye(2*Fem.NNode,2*Fem.NNode);
        KTtmp(pDof,:) = 0*I(pDof,:);
        KTtmp(pDof,pDof) = -EMod*eye(length(pDof));
        Ktr = KTtmp;

        if Fem.Iteration == 1
        F(pDof) = Fe(pDof)-EMod*F(pDof)*Fem.TimeStep*Fem.OptFactor;
        else, F(pDof) = Fe(pDof); 
        end
    else
        Ktmp = full(K); I = eye(2*Fem.NNode,2*Fem.NNode);
        Ktmp(pDof,:) = I(pDof,:);
        Ktmp(:,pDof) = I(:,pDof);
        F = -beta*K*F;
        F(pDof) = beta*Fem.Load(1:NLoad,id+1);
        K = Ktmp;
    end
end

if ~isempty(Fem.Contraction)
    Fb = sparse(Fem.i,Fem.e,Fem.fb0);
    F(:,1) = beta*Fb*Fem.Contraction(:);
end

if Fem.VolumetricPressure
    A = Fem.Mesh.get('Area');
    CellVolume = zeros(Fem.NElem,1);
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
    CellVolume(id) = A(id)*mean(Fem.PressureCell(:,2));
    Fb = sparse(Fem.i,Fem.e,Fem.fb0);
    W = sparse(Fem.i,Fem.j,Fem.ed);
    df0 = beta*W*Fb*CellVolume(:);
    F(:,1) = F(:,1) + df0;
    
    if Fem.SolverStartMMA
        z0 = Fem.Density;
        dz = 1e-2;
        Fem.dFdE = zeros(2*Fem.NNode,Fem.NElem);
        Pc = Fem.Mesh.get('Center');
        f0 = beta*sparse(Fem.i,Fem.e,Fem.fb0)*CellVolume(:);
        bnd = boundary(Pc(id,1),Pc(id,2));
        idmax = id(bnd);
        for ii = idmax
            ze = z0; ze(ii) = ze(ii) + dz;
            Fem.Density = ze;
            %[E,~,~,~] = MaterialField(Fem);
            id = FindElements(Fem,'FloodFill',Fem,ze);
            Fb = beta*sparse(Fem.i,Fem.e,Fem.fb0);
            CV = zeros(Fem.NElem,1);
            CV(id) = A(id)*mean(Fem.PressureCell(:,2));
            fe = Fb*CV(:);
            Fem.dFdE(:,ii) = ((fe - f0)/(dz));
        end
        
        Fem.Density = z0;
    end
end

Fem.fInternal = Fe; 
Fem.fExternal = F - spMat*Fem.Utmp;
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
function [Fe,Qe,Me,Ce,Ke,Kte,Svme,SS] = ...
    Locals(Fem,eNode,eDof)

nn = length(eNode);
Fe = zeros(2*nn,1);
Me = zeros(2*nn,2*nn);
Ce = zeros(2*nn,2*nn);
Ke = zeros(2*nn,2*nn);
Kte = zeros(2*nn,2*nn);
Qe = ones(nn,1);
W = Fem.ShapeFnc{nn}.W;
Q = Fem.ShapeFnc{nn}.Q;

SGP = zeros(length(W),6);
Nshp = length(Fem.ShapeFnc{nn}.N(:,:,1));
NNe = zeros(length(W),Nshp);

% quadrature loop
for q = 1:length(W)
    dNdxi = Fem.ShapeFnc{nn}.dNdxi(:,:,q);
    N = Fem.ShapeFnc{nn}.N(:,:,q);
    J0 = Fem.Node0(eNode,:).'*dNdxi;
    %Xe = Fem.Node(eNode,:);
    dNdx = dNdxi/J0;
    dJ = abs(det(J0));
       
    % get displacement field
    Delta = Fem.Utmp(eDof,:);
    
    % deformation gradient   
    F = DeformationGradient(Fem,Delta,dNdx);
    
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
    
    % linear strain-displacement operator
    [B,~,~] = LinearStrainOperator(Fem,N,dNdx);
    
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
    
    % true stress
    Scauchy = (1/det(F))*F*S0*(F.');
    
    SGP(q,:) = VoightNotation(Scauchy);
    
    % construct shaping
    NNe(((q-1)*Nshp + 1):(q*Nshp)) = N(:).';
end

SS = NNe.'*SGP;
[Svm, ~] = VonMises(SS(:,1),SS(:,2),SS(:,3),SS(:,4),SS(:,5),SS(:,6));
Svme = Svm(:); 
end

%----------------------------------------------------- deformation gradient
function F = DeformationGradient(Fem,U,dNdx)
nn = length(U)/2;
F = zeros(3);
UU = zeros(nn,2);
UU(:,1) = U(1:2:2*nn);
UU(:,2) = U(2:2:2*nn);
%F(1:2,1:2) = UU.'*dNdx;
F(1:2,1:2) = (dNdx'*UU)';
F = F + eye(3);
%if Fem.Nonlinear, F = F + eye(3); else, F = eye(3); end
end

%----------------------------------------------- TRIANGULAR SHAPE FUNCTIONS
function [Node,Node0] = UpdateNode(Fem,U)
Node0 = Fem.Node0; Ntmp = Node0;

[ux,uy,~] = DisplacementField(Fem,U);

Ntmp(:,1) = Node0(:,1) + ux(:);
Ntmp(:,2) = Node0(:,2) + uy(:);

Node = Ntmp;
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

%---------------------------------------------- COMPUTE CENTROID OF POLYGON
function [Pc,A] = ComputeCentroid(Fem)

Pc = zeros(Fem.NElem,2);
A  = zeros(Fem.NElem,1);

for el = 1:Fem.NElem
    
    % get nodal positions
    vx = Fem.Node(Fem.Element{el},1);
    vy = Fem.Node(Fem.Element{el},2);
    nv = length(Fem.Element{el});
    
    % shift the vertices by (+1)
    vxS = vx([2:nv 1]); vyS = vy([2:nv 1]);
    
    % compute volume of area
    tmp = vx.*vyS - vy.*vxS;
    A(el) = 0.5*sum(tmp);
    Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*tmp),...
        sum((vy+vyS).*tmp)];
end

end

%---------------------------------------------- COMPUTE CENTROID OF POLYGON
function n = ComputeNormal(Fem)

n{Fem.NElem} = [];

for ii = 1:Fem.NElem
    id = Fem.Element{ii};
    Nd = Fem.Node(id,:) - repmat(Fem.Center(ii,:),length(id),1);
    n{ii} = Nd;
end

end

%---------------------------------------------------------- strain operator
function [Bn,Bg,NN,tau] = NonlinearStrainOperator(Fem,N,dNdx,F)
nn = length(N);

NN = zeros(2,2*nn);
NN(1,1:2:2*nn) = N(:)';
NN(2,2:2:2*nn) = N(:)';
dNdxX = dNdx(:,1).';
dNdxY = dNdx(:,2).';

if strcmp(Fem.Type,'PlaneStress') || strcmp(Fem.Type,'PlaneStrain')
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
end

end 

%--------------------------------------------------- linear strain operator
function [B,NN,tau] = LinearStrainOperator(Fem,N,dNdx)

nn = length(N);

NN = zeros(2,2*nn); %XX = NN;
NN(1,1:2:2*nn) = N(:)';
NN(2,2:2:2*nn) = N(:)';
% XX(1,1:2:2*nn) = Xe(:,1);
% XX(1,2:2:2*nn) = Xe(:,2);

% X = XX*NN.';

if strcmp(Fem.Type,'PlaneStress') || strcmp(Fem.Type,'PlaneStrain')
    B = zeros(3,2*nn);
    B(1,1:2:2*nn) = dNdx(:,1)';
    B(2,2:2:2*nn) = dNdx(:,2)';
    B(3,1:2:2*nn) = dNdx(:,2)';
    B(3,2:2:2*nn) = dNdx(:,1)';
    tau = 1;
elseif strcmp(Fem.Type,'AxiSymmetric')
%     B = zeros(4,2*nn);
%     B(1,1:2:2*nn) = dNdx(:,1)';
%     B(2,2:2:2*nn) = dNdx(:,2)';
%     B(3,1:2:2*nn) = N.'/X(1);
%     B(4,1:2:2*nn) = dNdx(:,2)';
%     B(4,2:2:2*nn) = dNdx(:,1)';
%     tau = 2*pi*X(1);
end

end

%------------------------------------------------------ isotropic reduction
function [S, D, G] = IsotropicReduction(Fem,D0,S0)
D = zeros(3,3); S = zeros(3,1);

if strcmp(Fem.Type,'PlaneStrain')
D(1:2,1:2) = D0(1:2,1:2); D(3,3) = D0(4,4);
SIG = [S0(1),S0(4);
       S0(4),S0(2)];
G = kron(eye(2),SIG);
% G(1:2,1:2) = SIG;
% G(3:4,3:4) = SIG;
S(1) = S0(1); 
S(2) = S0(2); 
S(3) = S0(4); 

elseif strcmp(Fem.Type,'PlaneStress')
%D0 = inv(D0); 
D(1:2,1:2) = D0(1:2,1:2); D(3,3) = D0(4,4);
%D = inv(D);
SIG = [S0(1),S0(4);
       S0(4),S0(2)];
G = kron(eye(2),SIG);
% G(1:2,1:2) = SIG;
% G(3:4,3:4) = SIG;
S(1) = S0(1); 
S(2) = S0(2); 
S(3) = S0(4); 

elseif strcmp(Fem.Type,'AxiSymmetric')

end
end

%--------------------------------------------------- elemental force vector
function Ve = ElemCellVolumeForce(Fem,eNode,ElemId)

nn = length(eNode);
Ve = zeros(2*nn,1);

Norm = Fem.Normal;

for ii = 1:nn
    a = Norm{ElemId}(ii,1);
    b = Norm{ElemId}(ii,2);
    Ve(2*ii - 1,1) = a/norm([a,b]);
    Ve(2*ii,1) = b/norm([a,b]);
end

end

%------------------------------------------------------ isotropic reduction
function [flag,Fem] = CheckConvergence(Fem,R,SingularKt)

R = full(R);

if ((norm(R) > Fem.ResidualNorm) && (Fem.Iteration <= Fem.MaxIteration)  ...
        && ~SingularKt )
    %ProcessMonitor(Fem);
    flag = 0;
else
    
    if (Fem.Iteration > Fem.MaxIteration) || SingularKt
        flag = 2;
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
    Fem.TimeStep = max(Fem.TimeStep/2,Fem.TimeStepMin);
    Fem.Time = Fem.Time - Fem.TimeStep;
    Fem.MaxIteration = Fem.MaxIteration + 15;
end

Fem.Iteration = 1;

if ~Fem.SolverStart, Fem.Increment = 1; Fem.SolverStart = true;
else, Fem.Increment = Fem.Increment + 1;
end

Fem.LoadingFactor = (Fem.Time/Fem.TimeEnd);

if Fem.SolverStartMMA, ProcessMonitor(Fem); end

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

%----------------------------------------------------------- node adjecency
function A = NodeAdjecency(face)
n = max(cellfun(@max,face));      
A = sparse(n,n);
null = sparse(n,n);

for ii = 1:length(face)
    poly = double(face{ii});
    [a,b] = PermutationSet(numel(poly));
    B = null;
    B(poly(a),poly(b)) = 1;
    A = A + B;
end

A = double(A>0);
end

%----------------------------------------------------------- permuation set
function [i,j] = PermutationSet(n)
i = transpose(kron(1:n,ones(1,n)));
j = transpose(kron(ones(1,n),1:n));
end

%%/////////////////////////////////////////////////// TOPOLOGY OPTIMIZATION

%------------------------------------------------------ objective function
function [f,dfdE,dfdV] = ObjectiveFunction(Fem)
u(:,1) = Fem.Utmp;
fDof = GetFreeDofs(Fem);

if strcmp(Fem.OptimizationProblem,'Compliance') &&  ~Fem.Nonlinear
f = Fem.fExternal.'*u;
temp = cumsum(-u(Fem.i).*Fem.k.*u(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliance') &&  Fem.Nonlinear
E = MaterialField(Fem);
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
f = Fem.fExternal.'*u;
lam = 0*u(:,1);
lam(fDof) = K(fDof,fDof)\Fem.fExternal(fDof);
temp = cumsum(-u(Fem.i).*Fem.k.*lam(Fem.j));
%temp = cumsum(-Fem.fi(Fem.i).*lam(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliant') && ~Fem.Nonlinear
E = MaterialField(Fem);
f = Fem.OutputVector.'*u(:,1); 
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
lam = 0*u(:,1);
lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
temp = cumsum(-u(Fem.i,1).*Fem.k.*lam(Fem.j,1));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
elseif strcmp(Fem.OptimizationProblem,'Compliant') && Fem.Nonlinear
E = MaterialField(Fem);
f = Fem.OutputVector.'*u(:,1); 
K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
lam = 0*u(:,1);
lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
temp = cumsum(-u(Fem.i).*Fem.k.*lam(Fem.j));
%temp = cumsum(-Fem.fi(Fem.i).*lam(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
end

if Fem.VolumetricPressure && Fem.Nonlinear
% K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
% lam = 0*u(:,1);
% lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
% K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.k);
% lam = 0*u(:,1);
% lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
dfdE = dfdE(:) + (lam(:)'*Fem.dFdE)';
% temp = cumsum(Fem.dFdE(Fem.e).*lam(Fem.j));
% temp = temp(cumsum(Fem.ElemNDof.^2));
% dfdE = dfdE + [temp(1);temp(2:end)-temp(1:end-1)];
% elseif isfield(fem,'dFdE') && Fem.Nonlinear
%     dfdE = dfdE(:) + (u(:,2)'*fem.dFdE)';
elseif Fem.VolumetricPressure && ~Fem.Nonlinear
%     
temp = cumsum(Fem.dFdE(Fem.e).*lam(Fem.j));
temp = temp(cumsum(Fem.ElemNDof.^2));
dfdE = dfdE + [temp(1);temp(2:end)-temp(1:end-1)];
%     
%dfdE = dfdE(:) + (lam(:)'*Fem.dFdE)';
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
   
N = length(Fem.Density);
M = 1;

iter = Fem.IterationMMA; 
if isempty(Fem.fnorm), Fem.fnorm = abs(norm(f)); end

xval = Fem.Density;
xmin = zeros(N,1);
xmax = ones(N,1);
%f0val = (f/norm(f));
%df0dx = (dfdz/norm(f));
%f0val = (f/Fem.fnorm);
%df0dx = (dfdz/Fem.fnorm);
f0val = f;
df0dx = dfdz;
df0dx2 = 0*df0dx;
fval = g;
dfdx = dgdz;
dfdx2 = dgdz*0;

A0 = 1;
A = 0;
C = 10000*ones(M,1);
D = 0;

if iter == 1, Fem.low = xmin; Fem.upp = xmax; end
if iter > 1, Fem.xold1 = Fem.xold1; else, Fem.xold1 = 0; end
if iter > 2, Fem.xold2 = Fem.xold2; else, Fem.xold2 = 0; end

[xmma,~,~,~,~,~,~,~,~,Fem.low,Fem.upp] = mmasub(M,N,iter,xval,xmin,xmax,...
    Fem.xold1,Fem.xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,...
    Fem.low,Fem.upp,A0,A,C,D);

%zNew = xmma;

alpha = max(0.9-iter/180,0);
zNew = (1-alpha)*xmma + alpha*xval;

if iter >= 1, Fem.xold1 = xval; end
if iter >= 2, Fem.xold2 = Fem.xold1; end

end

%---------------------------------------------- method of moving asymptotes
function [Fem,zNew] = UpdateSchemeOC(Fem,dfdz,g,dgdz)  
move=0.1*(Fem.zMax-Fem.zMin); 
eta=0.5;
l1=0; l2=1e6;  
iter = Fem.IterationMMA; 
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

alpha = max(0.9-iter/80,0);
zNew = (1-alpha)*zNew + alpha*z0;

%Change = max(abs(zNew-z0))/(zMax-zMin);
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
function [E,dEdy,V,dVdy] = MaterialField(Fem)
y = Fem.SpatialFilter*Fem.Density;
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
d = cell(size(PS1,1),1);

for el = 1:size(PS1,1)   
    dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    
    if ~isempty(Fem.Periodic)
        dist2 = sqrt((PS1(el,1)-PS2(:,1) + Fem.Periodic(1)).^2 +...
            (PS1(el,2)-PS2(:,2) + Fem.Periodic(2)).^2);
        
        dist3 = sqrt((PS1(el,1)-PS2(:,1) - Fem.Periodic(1)).^2 +...
            (PS1(el,2)-PS2(:,2) - Fem.Periodic(2)).^2);
        
        dist = (min([dist,dist2,dist3].'))';
    end

    [I,J] = find(dist<=R);   
    d{el} = [I,J+(el-1),dist(I)];
end
% matrix of indices and distance value
d = cell2mat(d); 
end

end
end

%---------------------------------------------------------- material filter 
function ProcessMonitor(Fem)

FreeDofs = GetFreeDofs(Fem);

if Fem.SolverStartMMA
if Fem.IterationMMA == 1  && Fem.Increment == 1 && Fem.Iteration == 1
fprintf(' Iter | Inc  | Residual  | Obj. fun. | Vol.  | Change | p    |\n');
fprintf('--------------------------------------------------------------\n');
end

fprintf(' %1.0f\t  | %1.0f\t | %1.3e | %1.3e | %0.3f | %0.3f  | %0.2f \n',...
    Fem.IterationMMA,Fem.Increment,norm(Fem.Residual(FreeDofs)),...
    abs(Fem.Obj),Fem.Con+1,norm(Fem.Change),Fem.Penal);

else
if Fem.Iteration == 1  && Fem.Increment == 1
fprintf(' Iter | Inc  | Residual  | Max. Svm  | Vol.  | Change | p    |\n');
fprintf('--------------------------------------------------------------\n');
end

fprintf(' %1.0f\t  | %1.0f\t | %1.3e | %1.3e | %0.3f | %0.3f  | %0.2f \n',...
    Fem.Iteration,Fem.Increment,norm(Fem.Residual(FreeDofs)),...
    max(Fem.s(:,1)),Fem.Con+1,norm(Fem.Change),Fem.Penal);   
    
end
    
    

end

%------------------------------------------------------ optimization solver
function showInformation(Fem)
    
fprintf('==============================================================\n');
fprintf('Number of elem: %1.2fk \n',Fem.NElem/1e3);
fprintf('Element type: P4-P8 \n');
fprintf('Filter radius: %1.2f \n', Fem.FilterRadius);
fprintf('Max. iterations: %1.0f \n', Fem.MaxIterationMMA);
fprintf('Optimization: %s \n', Fem.OptimizationProblem);
fprintf('Material scheme: %s \n', Fem.MaterialInterpolation);
if Fem.Nonlinear, fprintf('Nonlinear geom: true\n');
else, fprintf('Nonlinear geom: false \n'); end
if strcmp(Fem.Material.Type,'Mooney'), fprintf('Material: Mooney \n');
c1 = Fem.Material.C10; c2 = Fem.Material.C01; d = Fem.Material.K;
fprintf('C10 = %1.e, C01 = %1.3e, K = %1.3e \n', c1,c2,d);
elseif strcmp(Fem.Material.Type,'Yeoh'), fprintf('Material = Yeoh \n');
c1 = Fem.Material.C1*1e3; 
c2 = Fem.Material.C2*1e3; 
c3 = Fem.Material.C3*1e3;
fprintf('C1 = %1.1e kPa, C2 = %1.1e kPa, C3 = %1.1e kPa\n', c1,c2,c3);
end

fprintf('==============================================================\n');
    
end

%--------------------------------------------------------------------------
function Rp = Reflection(Fem,P)
BdBox = Fem.BdBox;

Domain = @(x) dRectangle(x,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
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
