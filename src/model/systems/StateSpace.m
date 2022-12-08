classdef StateSpace

    properties (Access = public)
        NDim;         % state dimensions
        NInput;       % input dimension
        Nonlinear;    % nonlinear (1) or linear ODE (0)
        X0;           % initial conditions
        A;            % state matrix 
        B;            % input matrix
        C;            % output matrix
        Flow;         % system flow
        Hessian;      % global Hessian matrix
    end
    
    properties (Access = private)

    end
   
%--------------------------------------------------------------------------
methods  
%----------------------------------------------- MODAL SHAPE RECONSTRUCTION
function obj = StateSpace(Input,X0,varargin) 
    
    if ~isa(Input,'cell') && ~isa(Input,'function_handle')
       error(['Input must be a cell; for linear ODE: {A,B,C} or ',...
           'for nonlinear ode: {@f,@dfdq} or {@f}']); 
    end
    
    if  isa(Input,'function_handle') 
        Input = {Input};
    end

    % construct linear ODE
    if isa(Input{1},'double')
        obj.Nonlinear = false;
        obj.NDim = numel(X0(:));
        obj.X0   = X0(:);
        obj.A    = Input{1};
        obj.B    = Input{2};

        obj.NInput = size(obj.B,2);
        
        if size(obj.A,1) ~= numel(obj.X0)
            error(['Dimension error, matrix A,B have different'...
                ' dimension than X0']); 
        end

    elseif isa(Input{1},'function_handle')
        obj.Nonlinear = true;
        obj.NDim = numel(X0(:));
        obj.X0   = X0(:);

        obj.Flow = Input{1};

        if numel(Input) > 1
            f = Input{1};
            obj.Hessian = Input{2};
        end
    end

    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end
%---------------------------------------------------------------------- get     
function varargout = get(StateSpace,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Model.(varargin{ii});
        end
    else
        varargout = StateSpace.(varargin);
    end
end       
%---------------------------------------------------------------------- set
function StateSpace = set(StateSpace,varargin)
    for ii = 1:2:length(varargin)
        StateSpace.(varargin{ii}) = varargin{ii+1};
    end
end
%-------------------------------------------------------------- return flow
function dx = flow(StateSpace,x,varargin)
   if ~isempty(varargin)
       u = varargin{1}(:);
       if numel(varargin)>1
           t = varargin{2};
       end
   else
       t = 0;
       u = 0;
   end

   if ~StateSpace.Nonlinear
       dx = StateSpace.A*x(:) + StateSpace.B*u(:);
   else
       dx = StateSpace.Flow(x(:),u(:),t);
   end
end
%----------------------------------------------------------- return hessian
function H = hessian(StateSpace,x,varargin)
    
   if ~isempty(varargin)
       u = varargin{1}(:);
       if numel(varargin)>1
           t = varargin{2};
       end
   else
       t = 0;
       u = 0;
   end 


   if ~StateSpace.Nonlinear
       H = StateSpace.A;
   else
       if isempty(StateSpace.Hessian)
           eps = 1e-6;
           H = zeros(StateSpace.NDim);
           f0 = flow(StateSpace,x,u,t);
           dx = eye(StateSpace.NDim)*eps;
           for ii = 1:StateSpace.NDim
              f1 = flow(StateSpace,x + dx(:,ii),u,t);
              H(:,ii) = (f1-f0)/eps;
           end
       end
   end
end

end
methods (Access = private)
%---------------------------------------------------------------------- set

end
end

