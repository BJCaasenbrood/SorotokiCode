classdef Model

    properties (Access = public)
        NDim;         % Total state dimension
        NIn;          % Total input dimension
        NSys;         % Number of subsystems
        NState;       % Number of states per system
        NInput;       % Number of inputs per system
        TimeStep;     % Solver timesteps
        TimeEnd;      % Finite time horizon
        Systems;      % Cell of subsystems
        Controller;   % Global system controller
        Update;       % Update externals (pre-call to controller)
        Log;          % Logger
        
        t;            % Time
        X;            % States
        X0;           % Initial conditions
        U;            % Inputs
        Y;            % Outputs
        K;            % Linear controller
    end
    
    properties (Access = private)
        Flow;         % system flow
        Hessian;      % global Hessian matrix
        
        Xtmp;         % Temporary states
        Utmp;         % Temporary inputs
        
        MaxIteration;
        ResidualNorm;
        
        ShowProcess; 
        
    end
   
%--------------------------------------------------------------------------
methods  
%----------------------------------------------- MODAL SHAPE RECONSTRUCTION
function obj = Model(Input,varargin) 
    
    obj.NSys = 1;
    obj.MaxIteration = 20;
    obj.ResidualNorm = 1e-6;
    obj.TimeEnd  = 10;
    obj.TimeStep = 1/60;
    obj.ShowProcess = true;
    obj.Log = struct;
    
    if isa(Input,'Shapes')
        obj.NDim = Input.NDim;        
        obj.NIn  = Input.NInput;        
        obj.Systems{1} = Input;
        obj.NState{1}  = Input.NDim;
        obj.NInput{1}  = Input.NInput;
    else
        obj.NDim = Input.NDim;        
        obj.NIn  = Input.NInput;        
        obj.Systems{1} = Input;
        obj.NState{1}  = Input.NDim;
        obj.NInput{1}  = Input.NInput;
    end

    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    obj = obj.rebuild();
    
end
%---------------------------------------------------------------------- get     
function varargout = get(Model,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Model.(varargin{ii});
        end
    else
        varargout = Model.(varargin);
    end
end       
%---------------------------------------------------------------------- set
function Model = set(Model,varargin)
    for ii = 1:2:length(varargin)
        Model.(varargin{ii}) = varargin{ii+1};
    end
end
%---------------------------------------------------------------------- set
function StateList = findState(Model,Id,varargin)
    List = 1:Model.NDim;
    indexList =  cumsum(vertcat(Model.NState{:}));
    index = indexList(Id);
    ind1 = index+varargin{1};
    StateList = List(ind1);
end
%---------------------------------------------------------------------- set
function Model = addSystem(Model,varargin)
    for ii = 1:numel(varargin)
        Model.Systems{end+1} = varargin{1};
    end
    
    Model = rebuild(Model);
end
%------------------------------------------------------------------ rebuild
function Model = rebuild(Model)
    nDim = 0; nIn = 0;
    Xini = []; Uini = [];
    for ii = 1:numel(Model.Systems)
        Model.NState{ii} = Model.Systems{ii}.NDim;
        Model.NInput{ii} = Model.Systems{ii}.NInput;
        Xini = [Xini;Model.Systems{ii}.X0(:)];
        Uini = [Uini;zeros(Model.NInput{ii},1)];
        nDim = nDim + Model.Systems{ii}.NDim;
        nIn = nIn + Model.Systems{ii}.NInput;
    end
    
    Model.NSys = numel(Model.Systems);
    Model.NDim = nDim;
    Model.NIn  = nIn;
    Model.U    = Uini;
    Model.Log.x = [];
    Model.Log.t = [];
    Model.Log.u = [];
    
    if norm(Model.X0) == 0
      Model.X0   = Xini; 
    end
    
    [~, Model] = FlowEval(Model,Xini,0);
    
end
%------------------------------------------------------------------ rebuild
function Model = simulate(Model)

% quick rebuild    
Model = rebuild(Model);   
    
h  = Model.TimeStep;
nd = Model.NDim;
Ts = 0:h:Model.TimeEnd;
z  = Model.X0;

if Model.ShowProcess
    tic;
    disp('----------------------------------');
    fprintf('* Computing SR dynamics... ');
    progress('start')
end


% initial logger
Model.X = z(:);

for ii = 1:length(Ts)
    
    % assign states to temp. variable
    dX  = 1e3;
    x_  = z;
    itr = 1;
    
    while (dX > Model.ResidualNorm && itr <= Model.MaxIteration)
        
        Model.t = Ts(ii);
        
        if ~isempty(Model.Controller)
            Model.U = Model.Controller(Model);
        end

        % compute flow field
        if itr == 1
            [f1, Model] = FlowEval(Model,z,Model.t);
            dR = -h*f1;
        else
            [f2, Model] = FlowEval(Model,x_,Model.t+h);
            dR = x_ - z - 0.5*h*(f1 + f2);
        end
        
        % compute hessian
        Hes = HessianEval(Model,x_);
        
%         if ~isempty(Model.Controller)
%            HesU = ControllerHessianEval(Model,x_);
%         end
%         
        % hessian update iteration
        A  = -(-(1/2)*h*(Hes) + eye(nd,nd));
        dx = A\dR;

        % update local state solution
        x_  = x_  + dx;
        itr = itr + 1;
        
        Model.X = x_;
        
        % compute convergence 
        dX = norm(dx(1:Model.NDim));
        
        %
        if Model.ShowProcess
            inc = round(100*(ii/(length(Ts)))-1);
            progress(inc,100);
        end
    end
    
    % update states
    z = x_;

    Model.Log.t = [Model.Log.t,Model.t];
    Model.Log.x = [Model.Log.x;x_.'];
    Model.Log.u = [Model.Log.u;Model.U.'];
    
    if ~isempty(Model.Update)
        Model = Model.Update(Model);
    end  
end

if Model.ShowProcess
    progress(100,100);
    progress('end');
    disp('----------------------------------');
    tt = toc;
    disp(['* Number of Dof.   = ', num2str(Model.NDim)]);
    disp(['* Computation time = ', num2str(tt,3), ' s']);
    disp(['* Computation freq = ', num2str(round(1/h)), ' Hz']);
    disp(['* Real-time ratio  = ', num2str((Ts(end))/(tt),4), '']);
    
    disp('----------------------------------');
end

end
end

methods (Access = private)
%---------------------------------------------------------------------- set
function [F, Model] = FlowEval(Model,X,t)
F = zeros(Model.NDim,1);
Ui = Model.U;

indexQ = 0; indexU = 0;
for ii = 1:Model.NSys
    NDofQ = Model.NState{ii};
    NDofU = Model.NInput{ii};
    ind1 = indexQ+(1:NDofQ);  % state indexing
    ind2 = indexU+(1:NDofU);  % input indexing
    
    if isa(Model.Systems{ii},'Shapes') || ...
            isa(Model.Systems{ii},'RigidBody') 
        [F(ind1),Model.Systems{ii}]  = ...
            Model.Systems{ii}.flow(X(ind1),Ui(ind2),t);
    else
        F(ind1) = Model.Systems{ii}.flow(X(ind1),Ui(ind2),t);
    end

    indexQ = indexQ + NDofQ;
    indexU = indexU + NDofU;
end

end
%---------------------------------------------------------------------- set
function H = HessianEval(Model,X)
ti = Model.t;
Ui = Model.U;

H = zeros(Model.NDim,Model.NDim);

indexQ = 0; indexU = 0;
for ii = 1:Model.NSys
    NDofQ = Model.NState{ii};
    NDofU = Model.NInput{ii};
    ind1 = indexQ+1:indexQ+NDofQ;  % state indexing
    ind2 = indexU+1:indexU+NDofU;  % input indexing

    H(ind1,ind1) = Model.Systems{ii}.hessian(X(ind1),Ui(ind2),ti);

    indexQ = indexQ + NDofQ;
    indexU = indexU + NDofU;
end

end
%---------------------------------------------------------------------- set
function H = ControllerHessianEval(Model,X)
eps = 1e-6;

H  = zeros(Model.NDim);
ti = Model.t;

[f0, Model] = FlowEval(Model,X,ti);
dx = eye(Model.NDim)*eps;

% for ii = 1:Model.NDim
%     
%     [~, Model]  = FlowEval(Model,X + dx(:,ii),ti);
%     Model.U = Model.Controller(Model);
%     [f1, Model] = FlowEval(Model,X,ti);
%     
%     H(:,ii) = (f1-f0)/eps;
% 
% end

end
end
end
